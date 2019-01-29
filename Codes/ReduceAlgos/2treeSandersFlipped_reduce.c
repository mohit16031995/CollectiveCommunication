#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define RUNS 100
#define max_processors 1000
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective

struct TreeNode* leftTreeNode[max_processors];
struct TreeNode* rightTreeNode[max_processors];


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

void flipTree(struct TreeNode* root) {
	if (root == NULL)
		return;
	struct TreeNode* temp;
	temp = root->left_child;
	root->left_child = root->right_child;
	root->right_child = temp;
	flipTree(root->left_child);
	flipTree(root->right_child);	
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
	int rank,p,rootId = 0, index,cdone=0;
	int vrank, vpeer,peer;
	long int count,i,j,SIZE,CSIZE,k,logical_chunk_no;
	char *ptr;

	int CHUNK;

	int recvToLeft = -1;
    	int recvToRight = -1;

	int leftChildren = 0;
	int rightChildren = 0;

 	int leftPeers[2];
	int rightPeers[2];


	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], outmsg[SIZE];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
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
		flipTree(root2); //for sake of uniformity		
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

		flipTree(root2);
		temp_node = root2;
		while (temp_node->left_child != NULL) {
			temp_node = temp_node->left_child;
		}
		temp_node->left_child = newNode(no_of_process, 1);
		temp_node->left_child->parent = temp_node;
		
		
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
	CHUNK = atoi(argv[2]);
	CSIZE = SIZE/CHUNK; 
	SIZE = CSIZE*CHUNK;

	
	int *msg1 = malloc(SIZE*2 * sizeof(int));
	int *selfmsg = malloc(SIZE * sizeof(int));
	int *msg2 = msg1+SIZE;
	int ready[CHUNK];

	// for  recvs 
	MPI_Status stt;
	MPI_Request *req1= calloc(CHUNK*2,sizeof(MPI_Request));
	MPI_Request *req2 = req1+CHUNK;
	//printf("I am rank %d, my req1 is at %p\n",rank,req1);
	for (i=0;i<CHUNK*2;i++)
		req1[i] = MPI_REQUEST_NULL;
	// for  send
	MPI_Status sstt[CHUNK];
	MPI_Request sreq[CHUNK];
	for (i=0;i<CHUNK;i++)
		sreq[i] = MPI_REQUEST_NULL;
	
	if (rank != 0) {
		leftChildren = 0;
		rightChildren = 0;
		recvToLeft = leftTreeNode[rank]->parent->process_id;
		recvToRight = rightTreeNode[rank]->parent->process_id;

		/////////////////////
		
		//count no of childs
            struct TreeNode* t1 = leftTreeNode[rank];
            if (t1->left_child != NULL) {
				leftChildren = 1;
				leftPeers[0] = t1->left_child->process_id;
            }
            if (t1->right_child != NULL) {
				leftChildren = 2;
				leftPeers[1] = t1->right_child->process_id;
            }
            t1 = rightTreeNode[rank];
            if (t1->left_child != NULL) {
				rightChildren = 1;
				rightPeers[0] = t1->left_child->process_id;
            }
            if (t1->right_child != NULL) {
				rightChildren = 2;
				rightPeers[1] = t1->right_child->process_id;
            }
	
	}
	else {
		leftChildren = 1;
		rightChildren = 1;
		leftPeers[0] = root->process_id;
		rightPeers[0] = root2->process_id;
	}
	double t1,t2,res;
	for (i=0;i<RUNS;i++)
	{
		for (int ll=0;ll<SIZE;ll++) {
			selfmsg[ll] = 1;
			//if(rank==0) msg[i] = selfmsg[i];	
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		cdone=0; count=0;
		for (k=0;k<CHUNK;k++)
			ready[k]=0;
		t1 = MPI_Wtime();

		// set up all recv from left and right tree
		if (leftChildren)				//if not leafleft setup all even recvs
			for (j=0;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,leftPeers[0],j,MPI_COMM_WORLD,&req1[j]);
				if (leftChildren==2)
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,leftPeers[1],j,MPI_COMM_WORLD,&req1[j+CHUNK]);
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from left tree child %d\n", rank,j,leftPeers[1]);
			}
		if (rightChildren)				//if not leafRight setup all odd recvs
			for (j=1;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,rightPeers[0],j,MPI_COMM_WORLD,&req1[j]);
				if (rightChildren==2)
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,rightPeers[1],j,MPI_COMM_WORLD,&req1[j+CHUNK]);
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from right tree child %d\n", rank,j,rightPeers[1]);
			}

		// setup all sends in left and right leaves
		if (!leftChildren || !rightChildren)
			for (j=0;j<CHUNK;j++) 
			{	
				if (!leftChildren && !(j%2)) 
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,recvToLeft,j,MPI_COMM_WORLD,&sreq[count++]);
				if (!rightChildren && (j%2))
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,recvToRight,j,MPI_COMM_WORLD,&sreq[count++]);
			}

		//MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		int sent_as_leaf = count;                 
////		printf ("%d leaf role completed\n",rank);
		// check which recv finish and setup send for them
		int chunks_to_recv = ((CHUNK+1)/2)*leftChildren + (CHUNK/2)*rightChildren;
		//chunks_to_recv = CHUNK;
		//printf ("%d, left=%d right=%d, rank %d\n",chunks_to_recv,leftChildren,rightChildren,rank);
		if(rank!=0) 
		{
			while(cdone < (chunks_to_recv)) 
			{		
				MPI_Waitany(CHUNK*2, req1, &index, &stt);
				//printf("ok waiting, rank %d",rank);
				if(index == MPI_UNDEFINED) 
				{
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				logical_chunk_no = (index>=CHUNK) ? (index-CHUNK) : index;
////				printf ("\trank %d received chunk %d (actual %d)\n",rank,logical_chunk_no,index);
				
				double t11 = MPI_Wtime();
				j=index*CSIZE;
				
					for (k=logical_chunk_no*CSIZE;k<(logical_chunk_no+1)*CSIZE;k++) 
					{
						selfmsg[k] = (msg1[j]+selfmsg[k]);
						j++;
					}
				double t12 = MPI_Wtime()-t11;
////				printf("Calculation time 1 %d %1.9f\n", CSIZE, t12);
				ready[logical_chunk_no]++;
				if ( (logical_chunk_no%2) ? ready[logical_chunk_no]==rightChildren : ready[logical_chunk_no]==leftChildren )
				{	//forward this chunk ;
////					printf("chunk no. %d is ready to forward from rank %d\n",logical_chunk_no,rank);
					
					if (!(logical_chunk_no%2))
						MPI_Isend(selfmsg+logical_chunk_no*CSIZE,CSIZE,MPI_INT,recvToLeft,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
					else
						MPI_Isend(selfmsg+logical_chunk_no*CSIZE,CSIZE,MPI_INT,recvToRight,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
				}
				cdone++;
				//printf("cdone = %d, rank %d",cdone,rank);
			}
////			printf("done forwarding %d\n",rank);
		}
		else
			while(cdone < (chunks_to_recv)) 
			{		
				MPI_Waitany(CHUNK*2, req1, &index, &stt);
				if(index == MPI_UNDEFINED) 
				{
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				logical_chunk_no = (index>=CHUNK) ? (index-CHUNK)*CSIZE : index*CSIZE;
				double t11 = MPI_Wtime();
				j=index*CSIZE;
				
					for (k=logical_chunk_no;k<logical_chunk_no+CSIZE;k++) 
					{
						selfmsg[k] = (msg1[j]+selfmsg[k]);		
						j++;
					}
			
				double t12 = MPI_Wtime() - t11;
//				printf("Calculation time 2 %d %1.9f\n", CSIZE, t12);	
				cdone++;
			}



		MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		
////		printf("My work is done here, rank %d\n",rank);
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		//		for (i=0; i<SIZE; i++)
		//			printf ("%d  ",selfmsg[i]);

		if(rank==0)
		{
			printf("Run %ld time %1.9lf\n", i+1,res);
		//	for (int ii=0; ii<SIZE; ii++)
		//		printf ("%d  ",selfmsg[ii]);
		//	printf("\n");
		}

	}

	MPI_Finalize();
} 
