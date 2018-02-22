#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <liblsb.h>

#include "fompi.h"
#define MAX_NODE_NAME 20
#define CABINETS_PER_ROW 10
#define NCOORDS 6
#define STRCOORDS 256
#define ITERATIONS 100
#define WARMUP 20
#define MPICHECK(X) {if (X!=MPI_SUCCESS) {printf("Error on "#X"\n"); MPI_Finalize(); exit(-1); }}



typedef uint32_t coord_t;

int get_coords(int rank, char * argv1, coord_t * coords){

    //char tmp[40];
    //sprintf(tmp, "name%i\0", rank);
    //FILE * f = fopen(tmp, "r");
    FILE * f = fopen(argv1, "r");
    if (f==NULL) { printf("Error opening %s\n", argv1); return 0; }
    fscanf(f, "c%i-%ic%is%in%i", &coords[0], &coords[1], &coords[2], &coords[3], &coords[4]); 
    fclose(f);
    coords[5] = rank;
    
    //DEBUG
    //printf("[Rank %i] cabinet: (%i, %i); chassis: %i; slot: %i; node: %i\n", rank, coords[0], coords[1], coords[2], coords[3], coords[4]);

    return 1;
}

const char * get_level_name(int lvl){
    const char * str[] = {"ROW", "CABINET", "CHASSIS", "SLOT", "NODE", "PROCESS"};
    return str[lvl];
}

int coords2str(coord_t * coords, char * str, int len){
    //char strcoords[256];
    int w=0;
    for (int i=0; i<NCOORDS; i++) 
        w += snprintf(str+w, len-w, " %i", coords[i]);
            
    return w;

}

void bench(foMPI_Win win, int peer, int msg_size){
    
    void * buff = malloc(msg_size);
    memset(buff, 0, msg_size);
    
    foMPI_Win_lock(MPI_LOCK_EXCLUSIVE, peer, 0, win);
    /*WARMUP*/
    for (int t=0; t<WARMUP; t++){
        foMPI_Get(buff, msg_size, MPI_BYTE, peer, 0, msg_size, MPI_BYTE, win);
        foMPI_Win_flush(peer, win);
    }

    /*MEASURE*/

    for (int t=0; t<ITERATIONS; t++){       
        LSB_Res();
        foMPI_Get(buff, msg_size, MPI_BYTE, peer, 0, msg_size, MPI_BYTE, win);
        foMPI_Win_flush(peer, win);
        LSB_Rec(msg_size);
    }
    foMPI_Win_unlock(peer, win);    

    for (int i=0; i<msg_size/sizeof(int); i++) { if(((int *) buff)[i]!=peer) printf("error! expected: %i; found: %i\n", peer, ((int*) buff)[i]); }

    free(buff);

}

int main(int argc, char * argv[]){
    int pname_len, rank, csize;
    char pname[MPI_MAX_PROCESSOR_NAME];
    char nodename[MAX_NODE_NAME];
    memset(pname, MPI_MAX_PROCESSOR_NAME, '\0');
    memset(nodename, MAX_NODE_NAME, '\0');

    foMPI_Init(&argc, &argv);
    MPI_Get_processor_name(pname, &pname_len);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &csize);
 
    
    srand(rank*1000000);

    coord_t coords[NCOORDS]; /* rank coords */
    coord_t * allcoords; /* all coords (only at rank 0) */

    
    LSB_Init("daint_bench", 0);
    allcoords = (coord_t *) malloc(sizeof(coord_t)*NCOORDS*csize);
    


    if (argc<=2) { printf("Need path and message size!!!\n"); exit(-1); }
    if (!get_coords(rank, argv[1], coords)) { printf("Unable to get coordinates\n"); exit(-1); } 

    /* Parse input: message sizes */
    int max_msg_size = 0;
    int msg_size_count = argc - 2;
    int * msg_size = (int *) malloc(sizeof(int)*(msg_size_count));
    if (msg_size==NULL) {printf("Error on malloc\n"); exit(-1);}
    for (int i=0; i<msg_size_count; i++){
        msg_size[i] = atoi(argv[i+2]);
        if (max_msg_size<msg_size[i]) max_msg_size=msg_size[i];
    }
    if (max_msg_size==0) { printf("Cannot send 0byte message\n"); exit(-1);}
       
    /* Gather all nodes info at rank 0 */
    MPI_Allgather(coords, NCOORDS, MPI_INT, allcoords, NCOORDS, MPI_INT, MPI_COMM_WORLD);

   

    /* Window setup */
    foMPI_Win win;
    void * win_mem = malloc(max_msg_size);
    MPICHECK(foMPI_Win_create(win_mem, max_msg_size, sizeof(uint8_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win));

    MPI_Barrier(MPI_COMM_WORLD);
    foMPI_Win_fence(0, win);
    
    for (int i=0; i<max_msg_size/(sizeof(int)); i++) { ((int *)win_mem)[i] = rank; }

    //printf("Allocating window of %i bytes (total exps: %i)\n", max_msg_size, argc-2);

   
    
    if (rank==0){    
        printf("Received coords: \n");
        for (int i=0; i<csize; i++){
            printf("\tRank %i: ", i);
            for (int j=0; j<NCOORDS; j++) printf("%i ", allcoords[i*NCOORDS + j]);
            printf("\n");
        }
    }

    int * peers = (int *) malloc(sizeof(int)*csize);
    
    char strcoords[STRCOORDS];
    LSB_Set_Rparam_int("INITIATOR", rank);
    //coords2str(&allcords[rank*NCOORDS], strcoords, STRCOORDS);
    //LSB_Set_Rparam_string("INITIATOR_COORDS", strcoords);

    for (int q=NCOORDS-1; q>=0; q--){
        int c, peer;
        int peercount=0;

        for (int i=0; i<csize; i++){
            for (c=0; c<q && allcoords[i*NCOORDS + c]==coords[c]; c++){;}
            //found = c>q && (c==NCOORDS || allcoords[i*NCOORDS + c] != coords[c]); 
            if (c==q && allcoords[i*NCOORDS + c] != coords[c] && i>rank){
                peers[peercount++] = i;
            }
        }

        
        //if (peercount==0) { printf("No node found at %s level\n", get_level_name(q)); continue; }
          
        if (peercount>0){

            peer = peers[rand() % peercount];
    
   
            /* debug */
            coords2str(&allcoords[peer*NCOORDS], strcoords, STRCOORDS);
            printf("[Rank %i] Level: %s; Selected peer: %i; coords:%s\n", rank, get_level_name(q), peer, strcoords);
    
            /* Perform test */
            LSB_Set_Rparam_string("LEVEL", get_level_name(q));
            LSB_Set_Rparam_int("TARGET", peer);
            for (int im=0; im<msg_size_count; im++) bench(win, peer, msg_size[im]); 
        }
        MPI_Barrier(MPI_COMM_WORLD);   
    }



    //MPI_Barrier(MPI_COMM_WORLD);
    foMPI_Win_fence(0, win);

    foMPI_Win_free(&win);
    
    MPI_Barrier(MPI_COMM_WORLD);

    
    free(allcoords);
    LSB_Finalize();
    
    free(msg_size);
    free(peers);
    foMPI_Finalize(); 

    return 1;
}

