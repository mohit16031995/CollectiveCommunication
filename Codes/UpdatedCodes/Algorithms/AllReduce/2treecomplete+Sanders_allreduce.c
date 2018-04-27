#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define RUNS 100
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective
int p = 0;

#define max_processors 1000
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective

struct TreeNode* leftTreeNode[max_processors];
struct TreeNode* rightTreeNode[max_processors];

int funcRP(int rank) {
	if (rank==0)
		return 0;
	return (p - (p-rank)/2)%p;
}

void printGivenLevel(struct TreeNode* root, int level);
int height(struct TreeNode* node);
void printLevelOrder(struct TreeNode* root);
struct TreeNode* newNode(int data, int tree);
struct TreeNode* constructleft(int start_id, int last_id);
struct TreeNode* constructTree(int no_of_pe, int start_id);
void traverse(struct TreeNode* root);
struct TreeNode* mirror(struct TreeNode* root, int no_of_pe);
struct TreeNode* addNode (int node_id, struct TreeNode* root, int tree);
void addParentColor(int node_id, int tree, int color);

struct TreeNode
{
    int process_id;
    struct TreeNode *left_child;
    int leftColor, rightColor;
    struct TreeNode* parent;
    struct TreeNode *right_child;
};

void printLevelOrder(struct TreeNode* root)
{
    int h = height(root);
    int i;
    for (i=1; i<=h; i++)
        printGivenLevel(root, i);
}

/* Print nodes at a given level */
void printGivenLevel(struct TreeNode* root, int level)
{
    if (root == NULL)
        return;
    if (level == 1)
        printf("%d ", root->process_id);
    else if (level > 1)
    {
        printGivenLevel(root->left_child, level-1);
        printGivenLevel(root->right_child, level-1);
    }
}

/* Compute the "height" of a tree -- the number of
    nodes along the longest path from the root node
    down to the farthest leaf node.*/
int height(struct TreeNode* node)
{
    if (node==NULL)
        return 0;
    else
    {
        /* compute the height of each subtree */
        int lheight = height(node->left_child);
        int rheight = height(node->right_child);

        /* use the larger one */
        if (lheight > rheight)
            return(lheight+1);
        else return(rheight+1);
    }
}

struct TreeNode* newNode(int data, int tree)
{
  struct TreeNode* node = (struct TreeNode*)malloc(sizeof(struct TreeNode));
  node->process_id = data;
  node->left_child = NULL;
  node->right_child = NULL;
  node->leftColor = -1;
  node->rightColor = -1;
  node->parent = NULL;
  if (tree == 0) {
    leftTreeNode[data] = node;
  }
  else {
    rightTreeNode[data] = node;
  }
  return(node);
}

struct TreeNode* constructCompleteBinary(int start_id, int last_id) {
//    	printf("construct Complete binary called start_id = %d, last_id = %d\n", start_id, last_id);
	if (start_id < last_id) {
   	int mid = ceil((start_id+last_id) / 2.0);
    	struct TreeNode* root = newNode(mid,0);
    	root->left_child = constructCompleteBinary(start_id, mid-1);
    	if (root->left_child != NULL)
        	root->left_child->parent = root;
    	root->right_child = constructCompleteBinary(mid+1, last_id);
    	if (root->right_child != NULL)
        	root->right_child->parent = root;

    	return root;	
    }
    else if (start_id == last_id) {
        return newNode(start_id,0);
    }
    else {
        return NULL;
    }
}
struct TreeNode* constructTree(int no_of_pe, int start_id) {
    //printf("ConstructTree called with no_of_pe = %d, start_id = %d\n", no_of_pe, start_id);
	if (no_of_pe <= 0) {
        return NULL;
    }
    if (no_of_pe == 1) {
        return newNode(start_id, 0);
    }
    int h = ceil((log10(no_of_pe+2.0)) / log10(2.0));
    int root_id = pow(2,h-1)-1 + start_id;
    struct TreeNode* root = newNode(root_id,0);
    root->left_child = constructCompleteBinary(start_id, root_id-1);
    if (root->left_child != NULL)
        root->left_child->parent = root;
    root->right_child = constructTree(start_id+no_of_pe-1-root_id, root_id+1);
    if (root->right_child != NULL)
        root->right_child->parent = root;
    return root;
}
void traverse(struct TreeNode* root) {
   if (root == NULL) {
    return;
   }
   printf("%d ", root->process_id);
   traverse(root->left_child);
   traverse(root->right_child);
}
struct TreeNode* mirror(struct TreeNode* root, int no_of_pe) {
    if (root == NULL) {
        return NULL;
    }
    struct TreeNode* temp = newNode(no_of_pe - root->process_id+1,1);
    temp->left_child = mirror(root->left_child, no_of_pe);
    if (temp->left_child != NULL)
        temp->left_child->parent = temp;
    temp->right_child = mirror(root->right_child, no_of_pe);
    if (temp->right_child != NULL)
        temp->right_child->parent = temp;

    return temp;
}
void addParentColor(int node_id, int tree, int color) {
    //printf("addColor : %d, %d %d\n", node_id, tree, color);
    struct TreeNode* node;
    if (tree == 0) {
        node = leftTreeNode[node_id];
    }
    else {
        node = rightTreeNode[node_id];
    }
    struct TreeNode* parentNode = node->parent;
    if (parentNode == NULL) {
        //printf("%d ka parent null hai\n", node_id);
        return;
    }
    if (parentNode->left_child == node && parentNode->leftColor == -1) {
  //      printf("left child hai, ");
        parentNode->leftColor = color;
        if (parentNode->right_child != NULL) {
            //printf("right child %d hai\n", parentNode->right_child->process_id);
            parentNode->rightColor = !(color);
            addParentColor(parentNode->right_child->process_id, !tree, color);
        }
    }
    else if (parentNode->right_child == node && parentNode->rightColor == -1){
        //printf("right child hai, ");
        parentNode->rightColor = color;
        if (parentNode->left_child != NULL) {
            //printf("left child %d hai\n", parentNode->left_child->process_id);
            parentNode->leftColor = !color;
            addParentColor(parentNode->left_child->process_id, !tree, color);
        }
    }
}

int main(int argc,char *argv[]){
	int rank, index, cdone=0, cdone2 = 0;
	long int count, i, j, k, SIZE, CSIZE, logical_chunk_no;
	long int CHUNK;
	char* ptr;

	int parentLeft = -1;
	int parentRight = -1;

	int leftChildren = 0;
	int rightChildren = 0;

	int leftPeers[2];
	int rightPeers[2];

	int parentLeft2 = -1;
	int parentRight2 = -1;

	int leftChildren2 = 0;
	int rightChildren2 = 0;

	int leftPeers2[2];
	int rightPeers2[2];



	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], selfmsg[SIZE];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//	righteers[0] = p-1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=3){	 
		if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
		exit(0);
	}	

	
	int total_process, no_of_process;
	total_process = p;
	no_of_process = total_process - 1;
	
    struct TreeNode *root, *root2, *top_node;
	if (no_of_process % 2 == 0) {
		int h = ceil((log10(no_of_process+2.0)) / log10(2.0));
		if (no_of_process == pow(2,h)-2) {
			root = constructCompleteBinary(1, no_of_process);
			root2 = mirror(root, no_of_process);
		}
		else {
			root = constructTree(no_of_process,1);
			root2 = mirror(root, no_of_process);
		}
	}
	else {
		int temp = no_of_process-1;
		int h = ceil((log10(temp+2.0)) / log10(2.0));
		if (temp == pow(2,h)-2) {
			root = constructCompleteBinary(1, temp);
			root2 = mirror(root, temp);
		}
		else {
			root = constructTree(temp,1);
			root2 = mirror(root, temp);
			
		}

		struct TreeNode* temp_node = root;
		while (temp_node->right_child != NULL) {
			temp_node = temp_node->right_child;
		}
		
		temp_node->right_child = newNode(no_of_process, 0);
		temp_node->right_child->parent = temp_node;
		
		temp_node = root2;
		while (temp_node->right_child != NULL) {
			temp_node = temp_node->right_child;
		}
		
		temp_node->right_child = newNode(no_of_process, 1);
		temp_node->right_child->parent = temp_node;
		
	}
	top_node = (struct TreeNode*)malloc(sizeof(struct TreeNode));
	top_node->process_id = 0;
	top_node->left_child = root;
	top_node->right_child = root2;
	top_node->parent = NULL;
	root2->parent = top_node;
	root->parent = top_node;
	top_node->leftColor = -1;
	top_node->rightColor = -1;

	SIZE = strtol(argv[1], &ptr, 10);
	long int len = SIZE / (sizeof(int));
	CHUNK = atoi(argv[2]);
	CSIZE = len/CHUNK; 
	SIZE = CSIZE*CHUNK;

	int *msg1 = malloc((SIZE+1) * sizeof(int));
	int *selfmsg = malloc((SIZE+1) * sizeof(int));
	int *msg2 = malloc((SIZE+1) * sizeof(int));
	int ready[CHUNK];
	int *Reducedmsg = malloc((SIZE+1) * sizeof(int));
	// for  recvs 
	MPI_Status stt;
	MPI_Request *req= calloc(CHUNK*3,sizeof(MPI_Request));
	MPI_Request *req2 = req+CHUNK;
	//printf("I am rank %d, my req1 is at %p\n",rank,req1);
	for (i=0;i<CHUNK*3;i++)
		req[i] = MPI_REQUEST_NULL;
	// for  send
	MPI_Status sstt[CHUNK*4];
	MPI_Request sreq[CHUNK*4];

	for (i=0;i<CHUNK*4;i++)
		sreq[i] = MPI_REQUEST_NULL;

	double t1,t2,res;


	if (rank != 0) 						//if not root
	{
		parentLeft = rank / 2;
		parentRight = funcRP(rank);
		if (2*rank < p) {
			leftChildren = 1;
			leftPeers[0] = (2*rank);
		}
		if ((2*rank)+1 < p) {
			leftChildren = 2;
			leftPeers[1] = (2*rank)+1;
		}
		if (2*rank - p > 0) {
			rightChildren = 1;
			rightPeers[0] = (2*rank - p);
		}
		if (2*rank-p-1 > 0) {
			rightChildren = 2;
			rightPeers[1] = (2*rank-p-1);
		}
	}
	else
	{
		leftChildren = 1;
		rightChildren = 1;
		leftPeers[0] = 1;
		rightPeers[0] = p-1;
	}
	
	if (rank != 0) {
		leftChildren2 = 0;
		rightChildren2 = 0;
		parentLeft2 = leftTreeNode[rank]->parent->process_id;
		parentRight2 = rightTreeNode[rank]->parent->process_id;

		/////////////////////
		
		//count no of childs
            struct TreeNode* t1 = leftTreeNode[rank];
            if (t1->left_child != NULL) {
				leftChildren2 = 1;
				leftPeers2[0] = t1->left_child->process_id;
            }
            if (t1->right_child != NULL) {
				leftChildren2 = 2;
				leftPeers2[1] = t1->right_child->process_id;
            }
            t1 = rightTreeNode[rank];
            if (t1->left_child != NULL) {
				rightChildren2 = 1;
				rightPeers2[0] = t1->left_child->process_id;
            }
            if (t1->right_child != NULL) {
				rightChildren2 = 2;
				rightPeers2[1] = t1->right_child->process_id;
            }
	
	}
	else {
		leftChildren2 = 1;
		rightChildren2 = 1;
		leftPeers2[0] = root->process_id;
		rightPeers2[0] = root2->process_id;
	}
	//	double timings[2][50][515];
	//printf("rank%d leftChild %d leftA %d leftB %d rightChild %d rightA %d rightB %d\n", rank, leftChildren2, leftPeers2[0], leftPeers2[1], rightChildren2, rightPeers2[0], rightPeers2[1]);
	for (int ll=0;ll<SIZE;ll++) {
			selfmsg[ll] = 1;
			//if(rank==0) msg[i] = selfmsg[i];	
		}
	for (i=0;i<RUNS;i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
	
		cdone=0; count=0, cdone2=0;
		for (k=0;k<CHUNK;k++)
			ready[k]=0;
		t1 = MPI_Wtime();

		// set up all recv from left and right tree
		if (leftChildren)				//if not leafleft setup all even recvs
			for (j=0;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,leftPeers[0],j,MPI_COMM_WORLD,&req[j]);
				if (leftChildren==2)
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,leftPeers[1],j,MPI_COMM_WORLD,&req[j+CHUNK]);
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from left tree child %d\n", rank,j,leftPeers[1]);
			}
		if (rightChildren)				//if not leafRight setup all odd recvs
			for (j=1;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,rightPeers[0],j,MPI_COMM_WORLD,&req[j]);
				if (rightChildren==2)
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,rightPeers[1],j,MPI_COMM_WORLD,&req[j+CHUNK]);
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from right tree child %d\n", rank,j,rightPeers[1]);
			}

		// setup all sends in left and right leaves
		if (!leftChildren || !rightChildren)
			for (j=0;j<CHUNK;j++) 
			{	
				if (!leftChildren && !(j%2)) 
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentLeft,j,MPI_COMM_WORLD,&sreq[count++]);
				if (!rightChildren && (j%2))
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentRight,j,MPI_COMM_WORLD,&sreq[count++]);
			}
		
		
		int sent_as_leaf = count;                 
		int chunks_to_recv = ((CHUNK+1)/2)*leftChildren + (CHUNK/2)*rightChildren;
///Bcast
		if(rank!=0) {// if not root setup all recvs
		   for(j=0;j<CHUNK;j++) {
			// left tree				
			if(!(j%2)) MPI_Irecv(Reducedmsg+j*CSIZE,CSIZE,MPI_INT,parentLeft2,CHUNK+j,MPI_COMM_WORLD,&req[j+CHUNK*2]);			
			// right tree
			else MPI_Irecv(Reducedmsg+j*CSIZE,CSIZE,MPI_INT,parentRight2,CHUNK+j,MPI_COMM_WORLD,&req[j+CHUNK*2]);
		   }
		} 				


		if(rank!=0) 
		{
			while((cdone < (chunks_to_recv)) || (cdone2 < CHUNK)) 
			{		
				MPI_Waitany(CHUNK*3, req, &index, &stt);
				//printf("ok waiting, rank %d",rank);
				if(index == MPI_UNDEFINED) 
				{
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				if (index < CHUNK*2) {
					logical_chunk_no = (index>=CHUNK) ? (index-CHUNK) : index;
					j=logical_chunk_no*CSIZE;
					if (index < CHUNK) {					
						for (k=logical_chunk_no*CSIZE;k<(logical_chunk_no+1)*CSIZE;k++) 
						{
							selfmsg[k] = (msg1[j]+selfmsg[k]);
							j++;
						}
					}
					else {
						for (k=logical_chunk_no*CSIZE;k<(logical_chunk_no+1)*CSIZE;k++) 
						{
							selfmsg[k] = (msg2[j]+selfmsg[k]);
							j++;
						}
					}
					ready[logical_chunk_no]++;
					if ( (logical_chunk_no%2) ? ready[logical_chunk_no]==rightChildren : ready[logical_chunk_no]==leftChildren )
					{	//forward this chunk ;
						//printf("chunk no. %d is ready to forward from rank %d\n",logical_chunk_no,rank);
					
						if (!(logical_chunk_no%2))
							MPI_Isend(selfmsg+logical_chunk_no*CSIZE,CSIZE,MPI_INT,parentLeft,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
						else
							MPI_Isend(selfmsg+logical_chunk_no*CSIZE,CSIZE,MPI_INT,parentRight,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
					}
					cdone++;
					//printf("cdone = %d, rank %d",cdone,rank);
				}
				else {
					index = index - (CHUNK*2);
					if(index%2) {  // right tree
						if(rightChildren2) {
							MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,rightPeers2[0],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
						}
						if(rightChildren2==2) {
							MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,rightPeers2[1],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
						}
					} else {
						if(leftChildren2) {
							MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,leftPeers2[0],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
						}
						if(leftChildren2==2) {
							MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,leftPeers2[1],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
						}
					}
					cdone2++;
				}
			}
////			printf("done forwarding %d\n",rank);
		}
		else
			while(cdone < (chunks_to_recv)) 
			{		
				MPI_Waitany(CHUNK*2, req, &index, &stt);
				if(index == MPI_UNDEFINED) 
				{
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				logical_chunk_no = (index>=CHUNK) ? (index-CHUNK)*CSIZE : index*CSIZE;		
				j=logical_chunk_no;
				if (index < CHUNK) {
					for (k=logical_chunk_no;k<logical_chunk_no+CSIZE;k++) 
					{
						selfmsg[k] = (msg1[j]+selfmsg[k]);		
						j++;
					}
				}
				else {
					for (k=logical_chunk_no;k<logical_chunk_no+CSIZE;k++) 
					{
						selfmsg[k] = (msg2[j]+selfmsg[k]);		
						j++;
					}
				}
				cdone++;
				
				j=logical_chunk_no/CSIZE;
				ready[j]++;

				if ( (j%2) ? ready[j]==rightChildren : ready[j]==leftChildren ) {
				if(!(j%2)) { 
					if (leftChildren2)
						MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,leftPeers2[0],CHUNK+j,MPI_COMM_WORLD,&sreq[count++]);
					if (leftChildren2 == 2)
						MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,leftPeers2[1],CHUNK+j,MPI_COMM_WORLD,&sreq[count++]);
				}			
				// right tree
				else { 
					if (rightChildren2)
						MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,rightPeers2[0],CHUNK+j,MPI_COMM_WORLD,&sreq[count++]);
					if (rightChildren2 == 2)
						MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,rightPeers2[1],CHUNK+j,MPI_COMM_WORLD,&sreq[count++]);
				}					
				}
			}


		//printf("rank%d count%ld cdone=%d cdone2=%d\n", rank, count, cdone, cdone2);
		MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		//aprintf("DDDDD%d\n", rank);
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		if(rank==0)
		{
			printf("Run %ld time %1.9lf\n", i+1,res);
		}
		//else {
		//	for (int loop = 0; loop < SIZE; loop++) {
		//		printf("rank = %d %d\n", rank, Reducedmsg[loop]);
		//	}
		//}
	}
	//printf("out of loop RAnk = %d\n", rank);
	MPI_Finalize();
} 


