#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define max_processors 1000

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
int main(int argc, char *argv[]) {
  	int seed = time(NULL);
	srand(seed);
	
	int rank,p,index,cdone=0;
	long int sent,i,j,SIZE,CSIZE;
	char *ptr;

	int CHUNK;
	MPI_Init(&argc,&argv);
	//printf("point0.2\n");
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double t01, t02, t03, t_construction, t_coloring, t_testing = 0;
	t01 = MPI_Wtime(); 
	
	 int total_process, no_of_process;
	 total_process = atoi(argv[1]);
	no_of_process = total_process - 1;
	
	double timings[2][50][515];
	// receive 0 and send 1
	memset(timings, 0, sizeof(timings));
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

	t02 = MPI_Wtime();
	/*if (rank == 0) {	
    		printLevelOrder(root);
    		printf("\n");
    		printLevelOrder(root2);
    		printf("\n");
		printLevelOrder(top_node);
		printf("\n");
	}*/
	// coloring

	if (top_node->leftColor == -1) {
		top_node->leftColor = 0;
		addParentColor(top_node->left_child->process_id, 1, 1);
	}
	if (top_node->rightColor == -1) {
		top_node->rightColor = 1;
		addParentColor(top_node->right_child->process_id, 0, 0);
	}
	for (int i = 1; i <= no_of_process; i++) {
		struct TreeNode* temp = leftTreeNode[i];
		if (temp->left_child != NULL && temp->leftColor == -1) {
			temp->leftColor = 0;
			addParentColor(temp->left_child->process_id, 1, 1);
		}
		if (temp->right_child != NULL && temp->rightColor == -1) {
			temp->rightColor = 1;
			addParentColor(temp->right_child->process_id, 0, 0);
		}
		struct TreeNode* temp2 = rightTreeNode[i];
		if (temp2->left_child != NULL && temp2->leftColor == -1) {
			temp2->leftColor = 0;
			addParentColor(temp2->left_child->process_id, 0, 1);
		}
		if (temp2->right_child != NULL && temp2->rightColor == -1) {
			temp2->rightColor = 1;
			addParentColor(temp2->right_child->process_id, 0, 0);
		}
	}

	t03 = MPI_Wtime();	
	t03 = t03 - t02;
	t02 = t02 - t01;

	
// RUN MPI
    
	if(argc!=4){
	  if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
		printf("error due to input\n");
	  exit(0);
	}

	SIZE = strtol(argv[2], &ptr, 10);
	CHUNK = atoi(argv[3]);
	CSIZE = SIZE/CHUNK;
	SIZE = CSIZE*CHUNK;

	char *msg = malloc(SIZE+1 * sizeof(char));
	char *outmsg = malloc(SIZE+1 * sizeof(char));

	// for recvs
	MPI_Status stt;
	MPI_Request req[CHUNK];

	// for send
	MPI_Status sstt[CHUNK*2];
	MPI_Request sreq[CHUNK*2];

	double t1,t2,res;

	for(i=0;i<SIZE;i++) {
		outmsg[i] = 'A'+i%26;
		if(rank==0) msg[i] = outmsg[i];
	}
	outmsg[SIZE] = '\0';
	msg[SIZE] = '\0';

	int no_of_childs = 0;
	int left_tree_childs = 0;
	int right_tree_childs = 0;
	int left_child_rank = -1, right_child_rank = -1;
	if (rank != 0) {
	//count no of childs
            struct TreeNode* t1 = leftTreeNode[rank];
            if (t1->left_child != NULL) {
                no_of_childs += 1;
		left_tree_childs += 1;
		left_child_rank = t1->left_child->process_id;
            }
            if (t1->right_child != NULL) {
                no_of_childs += 1;
		left_tree_childs += 1;
		right_child_rank = t1->right_child->process_id;
		
            }
            t1 = rightTreeNode[rank];
            if (t1->left_child != NULL) {
                no_of_childs += 1;
		right_tree_childs += 1;
		left_child_rank = t1->left_child->process_id;
            }
            if (t1->right_child != NULL) {
                no_of_childs += 1;
		right_tree_childs += 1;
		right_child_rank = t1->right_child->process_id;
            }
	}
	
	MPI_Reduce(&t02, &t_construction, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&t03, &t_coloring, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	if (rank == 0) {
		printf("Tree Construction time = %1.9lf\n", t_construction);
		printf("Tree Coloring time = %1.9lf\n", t_coloring);
	}

    int RUNS = 1;
    for(i=0;i<RUNS;i++){
		MPI_Barrier(MPI_COMM_WORLD);
		
		t1 = MPI_Wtime();

        sent = 0;
        if (rank == 0){  // if root then setup all send's
            for(j=0;j<CHUNK;j++) {
			// left tree
                if(!(j%2)) {
                    MPI_Isend(outmsg+j*CSIZE,CSIZE,MPI_CHAR, root->process_id,j,MPI_COMM_WORLD, &sreq[sent++]);
			timings[1][j][root->process_id] = MPI_Wtime()-t1;
                }
			// right tree
                else {
                    MPI_Isend(outmsg+j*CSIZE,CSIZE,MPI_CHAR,root2->process_id,j,MPI_COMM_WORLD, &sreq[sent++]);
			timings[1][j][root2->process_id] = MPI_Wtime()-t1;
                }
            }
	}
        if(rank!=0) {// if not root setup all recvs
	
		MPI_Irecv(msg, CSIZE, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &req[0]);
		MPI_Irecv(msg+CSIZE, CSIZE, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &req[1]);
  //    	printf("Wait: rank = %d\n", rank);
		MPI_Waitany(2, req, &index, &stt);
            if(index == MPI_UNDEFINED) {
                printf("Unexpected error!\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
	//	printf("Recieved: Process id %d wait over recv %d chunk\n", rank, stt.MPI_TAG);
            int turn;
            int leftrecvd = -2;
            int rightrecvd = -1;
		int leftleftsent = 0;
		int leftrightsent = 0;
		int rightleftsent = 0;
		int rightrightsent = 0;
            if (index == 0) {
                leftrecvd = 0;
		timings[0][0][leftTreeNode[rank]->parent->process_id] = MPI_Wtime()-t1;
                struct TreeNode* temp = leftTreeNode[rank];
                if (temp->parent->left_child == temp) {
                    turn = temp->parent->leftColor;
                } 
                else {
                    turn = temp->parent->rightColor;
                }
            }
            else if (index == 1) {
                rightrecvd = 1;
		timings[0][1][rightTreeNode[rank]->parent->process_id] = MPI_Wtime()-t1;
                struct TreeNode* temp = rightTreeNode[rank];
                if (temp->parent->left_child == temp) {
                    turn = temp->parent->leftColor;
                }
                else {
                    turn = temp->parent->rightColor;
                }
            }
            else {
                printf("error\n");
            }
            int recieved = 1;
            
		
            
	
            int checkflag = 1;
	    int oddChunks = CHUNK / 2;
	    int evenChunks = CHUNK - oddChunks;
            while(recieved < CHUNK || sent < (left_tree_childs*evenChunks + right_tree_childs*oddChunks )) {
//		printf("Process %d entered loop, turn = %d\n", rank, turn);
		checkflag = 0;
                if (turn == 0) {
			turn = 1;
		}
		else {
			turn = 0;
		}

                if (leftrecvd != -2) {// left send
                    struct TreeNode* temp = leftTreeNode[rank];
                    if(temp->left_child != NULL && temp->leftColor == turn && leftleftsent <= leftrecvd/2) {
                        MPI_Isend(msg+leftrecvd*CSIZE,CSIZE,MPI_CHAR,temp->left_child->process_id,leftrecvd,MPI_COMM_WORLD,&sreq[sent++]);
			leftleftsent++;
			checkflag = 1;
	//		printf("%d rank process sent %d chunk to %d rank process-a\n", rank, leftrecvd, temp->left_child->process_id);
			timings[1][leftrecvd][temp->left_child->process_id] = MPI_Wtime()-t1;
                    }
                    else if (temp->right_child != NULL && temp->rightColor == turn && leftrightsent <= leftrecvd/2) {
                        MPI_Isend(msg+leftrecvd*CSIZE,CSIZE,MPI_CHAR,temp->right_child->process_id,leftrecvd,MPI_COMM_WORLD,&sreq[sent++]);
			leftrightsent++;
			checkflag =1;
			timings[1][leftrecvd][temp->right_child->process_id] = MPI_Wtime()-t1;
	//		printf("%d rank process sent %d chunk to %d rank process-b\n", rank, leftrecvd, temp->right_child->process_id);
                    }
                }
                if (rightrecvd != -1) { // right send
                    struct TreeNode* temp = rightTreeNode[rank];
                    if(temp->left_child != NULL && temp->leftColor == turn && rightleftsent <= rightrecvd/2) {
                        MPI_Isend(msg+rightrecvd*CSIZE,CSIZE,MPI_CHAR,temp->left_child->process_id,rightrecvd,MPI_COMM_WORLD,&sreq[sent++]);
			rightleftsent++;
			checkflag = 1;
			timings[1][rightrecvd][temp->left_child->process_id] = MPI_Wtime()-t1;
	//		printf("%d rank process sent %d chunk to %d rank process-c\n", rank, rightrecvd, temp->left_child->process_id);
                    }
                    else if (temp->right_child != NULL && temp->rightColor == turn && rightrightsent <= rightrecvd/2) {
                        MPI_Isend(msg+rightrecvd*CSIZE,CSIZE,MPI_CHAR,temp->right_child->process_id,rightrecvd,MPI_COMM_WORLD,&sreq[sent++]);
			rightrightsent++;
			checkflag = 1;
			timings[1][rightrecvd][temp->right_child->process_id] = MPI_Wtime()-t1;
//			printf("%d rank process sent %d chunk to %d rank process-d\n", rank, rightrecvd, temp->right_child->process_id);
                    }
                }

                // left tree
                if ((leftTreeNode[rank]->parent->left_child != NULL && leftTreeNode[rank]->parent->left_child == leftTreeNode[rank] && leftTreeNode[rank]->parent->leftColor == turn)
                 || (leftTreeNode[rank]->parent->right_child != NULL && leftTreeNode[rank]->parent->right_child == leftTreeNode[rank] && leftTreeNode[rank]->parent->rightColor == turn)) {
                   if (leftrecvd != -2) {
			 j = leftrecvd + 2;
                    if (j < CHUNK) {
	//		printf("Waiting: process %d waiting for packet %ld from %d\n", rank, j, leftTreeNode[rank]->parent->process_id);
                        int ierr = MPI_Recv(msg+j*CSIZE,CSIZE,MPI_CHAR,leftTreeNode[rank]->parent->process_id,j,MPI_COMM_WORLD,&stt);
                        if (ierr == MPI_SUCCESS) {
	//		    printf("Recieved: Process %d recieved packet %ld from %d\n", rank, j, leftTreeNode[rank]->parent->process_id);	
			    recieved += 1;
                            leftrecvd = j;
			    timings[0][j][leftTreeNode[rank]->parent->process_id] = MPI_Wtime()-t1;
				
                        }
                        else {
                            printf("error recieving chunk no. %ld from left tree with node %d\n", j, rank);
                        }
			checkflag = 1;
                    }
			}
		   else {
			int flag;
			MPI_Test(&req[0], &flag, &stt);
			if (flag) {
				recieved += 1;
				leftrecvd = 0;
			    timings[0][0][leftTreeNode[rank]->parent->process_id] = MPI_Wtime()-t1;
	//			printf("%d rank process recvd 0 chunk", rank);
			}	
			checkflag = 1;
	       	   }
                }
                if ((rightTreeNode[rank]->parent->left_child != NULL && rightTreeNode[rank]->parent->left_child == rightTreeNode[rank] && rightTreeNode[rank]->parent->leftColor == turn)
                 || (rightTreeNode[rank]->parent->right_child != NULL && rightTreeNode[rank]->parent->right_child == rightTreeNode[rank] && rightTreeNode[rank]->parent->rightColor == turn)) {
                    if (rightrecvd != -1) {
			j = rightrecvd + 2;
                    if (j < CHUNK) {
	//		printf("Waiting: process %d waiting for packet %ld from %d\n", rank, j, rightTreeNode[rank]->parent->process_id);
                       int ierr =  MPI_Recv(msg+j*CSIZE,CSIZE,MPI_CHAR,rightTreeNode[rank]->parent->process_id,j,MPI_COMM_WORLD,&stt);
                        if (ierr == MPI_SUCCESS) {
          //                  printf("Recieved: Process %d recved packet %ld from %d\n", rank, j, rightTreeNode[rank]->parent->process_id);
				recieved += 1;
                            rightrecvd = j;
			    timings[0][j][rightTreeNode[rank]->parent->process_id] = MPI_Wtime()-t1;
                        }
                        else {
                            printf("error recieving chunk no. %ld from right tree with node %d\n", j, rank);
                        }
			checkflag = 1;
                    }
		   }
		   else {
			double t_temp = MPI_Wtime();
			int flag;
			MPI_Test(&req[1], &flag, &stt);
			if (flag) {
				recieved += 1;
				rightrecvd = 1;
			    	timings[0][1][rightTreeNode[rank]->parent->process_id] = MPI_Wtime()-t1;
	//			printf("%d rank process recvd 1 chunk\n", rank);
			}	
			checkflag = 1;
			t_temp = t_temp - MPI_Wtime();
			t_testing += t_temp; 
		   }
                }
            } 
        }      
        MPI_Waitall(sent,sreq,sstt);  // wait for all send to finish
	//printf("Computation\n");
		double total_test;
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&t_testing, &total_test, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
		if(rank==0){
		
			printf("Run %ld time %1.9lf\n", i+1,res);
			printf("Process %d Run %ld testing time %1.9lf\n", rank, i+1,total_test);
			
			
		} else {
	//		printf("Outmsg %s Inmsg %s\n", outmsg,msg);
			j=strcmp(outmsg,msg);
			j!=0?printf("Error - msgs different\n"):0;
			memset(msg,'$',SIZE);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != 0) {
		for (int ll = 0; ll < CHUNK; ll++) {
			int from;
			if (ll % 2 == 0) {
				from = leftTreeNode[rank]->parent->process_id;
			}
			else {
				from = rightTreeNode[rank]->parent->process_id;
			}
			if (left_child_rank != -1 && right_child_rank != -1 && timings[1][ll][left_child_rank] != 0.0 && timings[1][ll][right_child_rank] != 0.0) {
				printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf, left_sent %d %1.9lf, right_sent %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from], left_child_rank, timings[1][ll][left_child_rank], right_child_rank, timings[1][ll][right_child_rank]);
			}
			else if(left_child_rank != -1 && timings[1][ll][left_child_rank] != 0.0) {
				printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf, left_sent %d %1.9lf\n", rank, i+1, ll, leftTreeNode[rank]->parent->process_id, timings[0][ll][leftTreeNode[rank]->parent->process_id], left_child_rank, timings[1][ll][left_child_rank]);
			}
			else if (right_child_rank != -1 && timings[1][ll][right_child_rank] != 0.0) {
				printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf, right_sent %d %1.9lf\n", rank, i+1, ll, leftTreeNode[rank]->parent->process_id, timings[0][ll][leftTreeNode[rank]->parent->process_id], right_child_rank, timings[1][ll][right_child_rank]);
			}
			else {
				printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from]);
			}
			//printf("\n");
		}
		}
		else {
			for(int ll=0;ll<CHUNK;ll++) {
			// left tree
				if(!(ll%2)) {
					printf("Logs, Process %d, Run %ld, chunk %d, left_sent %d %1.9lf\n", rank, i+1, ll, root->process_id, timings[1][ll][root->process_id]);
				}
					// right tree
				else {
				   printf("Logs, Process %d, Run %ld, chunk %d, right_sent %d %1.9lf\n", rank, i+1, ll, root2->process_id, timings[1][ll][root2->process_id]);
				}
            		}
		}
    } // End of loop for runs
    MPI_Finalize();
	
    return 0;
}
