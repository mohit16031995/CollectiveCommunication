#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define max_processors 1000

struct TreeNode* leftTreeNode[max_processors];
struct TreeNode* rightTreeNode[max_processors];

struct TreeNode
{
    int process_id;
    struct TreeNode *left_child;
    int leftColor, rightColor;
    struct TreeNode* parent;
    struct TreeNode *right_child;
};

void printGivenLevel(struct TreeNode* root, int level);
int height(struct TreeNode* node);

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
  if (tree == 0) {
    leftTreeNode[data] = node;
  }
  else {
    rightTreeNode[data] = node;
  }
  return(node);
}

struct TreeNode* constructleft(int start_id, int last_id) {
    if (start_id < last_id) {
    int mid = (last_id+start_id) / 2;
    struct TreeNode* root = newNode(mid,0);
    root->left_child = constructleft(start_id, mid-1);
    if (root->left_child != NULL)
        root->left_child->parent = root;
    root->right_child = constructleft(mid+1, last_id);
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
    //printf("construct %d %d \n" , no_of_pe, start_id);
    if (no_of_pe <= 0) {
        return NULL;
    }
    if (no_of_pe == 1) {
        return newNode(start_id,0);
    }
    int h = ceil((log10(no_of_pe+2.0)) / log10(2.0));
    int root_id = pow(2,h-1)-1 + start_id;
    struct TreeNode* root = newNode(root_id,0);
    //printf("%d", root->process_id);
    root->left_child = constructleft(start_id, root_id-1);
    if (root->left_child != NULL)
        root->left_child->parent = root;
    root->right_child = constructTree(no_of_pe-root_id-1, root_id+1);
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
    struct TreeNode* temp = newNode(no_of_pe - 1- root->process_id,1);
    temp->left_child = mirror(root->left_child, no_of_pe);
    if (temp->left_child != NULL)
        temp->left_child->parent = temp;
    temp->right_child = mirror(root->right_child, no_of_pe);
    if (temp->right_child != NULL)
        temp->right_child->parent = temp;

    return temp;
}
struct TreeNode* addNode (int node_id, struct TreeNode* root, int tree) {
    if (root == NULL) {
        return newNode(node_id, tree);
    }
    else {
        root->right_child = addNode(node_id, root->right_child, tree);
        return root;
    }
}
void addParentColor(int node_id, int tree, int color) {
    printf("addColor : %d, %d %d\n", node_id, tree, color);
    struct TreeNode* node;
    if (tree == 0) {
        node = leftTreeNode[node_id];
    }
    else {
        node = rightTreeNode[node_id];
    }
    struct TreeNode* parentNode = node->parent;
    if (parentNode == NULL) {
        printf("%d ka parent null hai\n", node_id);
        return;
    }
    if (parentNode->left_child == node) {
        printf("left child hai, ");
        parentNode->leftColor = color;
        if (parentNode->right_child != NULL) {
            printf("right child %d hai\n", parentNode->right_child->process_id);
            parentNode->rightColor = !(color);
            addParentColor(parentNode->right_child->process_id, !tree, color);
        }
    }
    else if (parentNode->right_child == node){
        printf("right child hai, ");
        parentNode->rightColor = color;
        if (parentNode->left_child != NULL) {
            printf("left child %d hai\n", parentNode->left_child->process_id);
            parentNode->leftColor = !color;
            addParentColor(parentNode->left_child->process_id, !tree, color);
        }
    }
}
int main() {
    int no_of_process;
    scanf("%d", &no_of_process);
    struct TreeNode* root, *root2;
    int h = ceil((log10(no_of_process)) / log10(2.0));
    if (no_of_process == pow(2,h) - 2) {
        root = constructleft(0, no_of_process-1);
        root2 = mirror(root, no_of_process);
    }
    else if (no_of_process == pow(2,h)-1) {
        root = constructleft(0, no_of_process-2);
        root2 = mirror(root, no_of_process-1);
        struct TreeNode* temp = newNode(no_of_process-1,0);
        temp->left_child = root;
        root = temp;
        struct TreeNode* temp2 = newNode(no_of_process-1,1);
        temp2->left_child = root2;
        root2 = temp2;
    }
    else if (no_of_process % 2 == 0) {
        root = constructTree(no_of_process, 0);
        root->parent = NULL;
        root2 = mirror(root, no_of_process);
        root2->parent = NULL;
    }
    else {
        root = constructTree(no_of_process-1, 0);
        root->parent = NULL;
        root2 = mirror(root, no_of_process-1);
        root2->parent = NULL;
        addNode(no_of_process-1, root, 0);
        addNode(no_of_process-1, root2, 1);

    }
    printf("Constructed \n");
    printLevelOrder(root);
    printf("\n");
    printLevelOrder(root2);
    printf("\n");

    if (root->leftColor == -1 && root->rightColor == -1) {
        if (root->left_child != NULL) {
            root->leftColor = 0;
            addParentColor(root->left_child->process_id, 1, 1);
        }
        if (root->right_child != NULL) {
            root->rightColor = 1;
            addParentColor(root->right_child->process_id, 1, 0);
        }
    }
    if (root2->leftColor == -1 && root2->rightColor == -1) {
        if (root2->left_child != NULL) {
            root2->leftColor = 0;
            addParentColor(root2->left_child->process_id, 0, 1);
        }
        if (root2->right_child != NULL) {printLevelOrder(root);
        printf("\n");
        printLevelOrder(root2);
        printf("\n");
            root2->rightColor = 1;
            addParentColor(root2->right_child->process_id, 0, 0);
        }
    }
    printf("Colored\n");
    for (int i = 0; i < no_of_process; i++) {
        printf("process_id = %d : \n", i);
        struct TreeNode* temp1 = leftTreeNode[i];
        struct TreeNode* temp2 = rightTreeNode[i];

        printf("Tree1 : leftColor = %d, rightColor = %d\n", temp1->leftColor, temp1->rightColor);
        printf("Tree2 : leftColor = %d, rightColor = %d\n\n", temp2->leftColor, temp2->rightColor);
    }

    // RUN MPI
    int rank,p,index,cdone=0;
    long int count,i,j,SIZE,CSIZE;
	char *ptr;

	int CHUNK;

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

	SIZE = strtol(argv[1], &ptr, 10);
	CHUNK = atoi(argv[2]);
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
		outmsg[i] = 'A'+rand()%26;
		if(rank==root) msg[i] = outmsg[i];
	}
	outmsg[SIZE] = '\0';
	msg[SIZE] = '\0';

    for(i=0;i<RUNS;i++){
		MPI_Barrier(MPI_COMM_WORLD);
        count = 0;
        if (rank == no_of_process){  // if root then setup all send's
		   for(j=0;j<CHUNK;j++) {
			// left tree
			if(!(j%2)) {
				MPI_Isend(outmsg+j*CSIZE,CSIZE,MPI_CHAR, root->process_id,j,MPI_COMM_WORLD,&sreq[count++]);
			}
			// right tree
			else {
				MPI_Isend(outmsg+j*CSIZE,CSIZE,MPI_CHAR,root2->process_id,j,MPI_COMM_WORLD,&sreq[count++]);
			}
		   }
		}
        if(rank!=no_of_process) {// if not root setup all recvs

          int recieved = 0;
          int sent = 0;
          bool turn = 1;
          int leftchunkrecvd = -2;
          int rightchunkrecvd = -1;
          
		  while(recieved < CHUNK || sent < CHUNK) {
            turn = turn + 1;

            if (leftchunkrecvd != -2) {// left send
                struct TreeNode* temp = leftTreeNode[rank];
                if(temp->left_child != NULL && temp->leftColor == turn) {
                    MPI_Isend(msg+leftchunkrecvd*CSIZE,CSIZE,MPI_CHAR,temp->left_child->process_id,leftchunkrecvd,MPI_COMM_WORLD,&sreq[count++]);
                    sent += 1;
                }
                else if (temp->right_child != NULL && temp->rightColor == turn) {
                    MPI_Isend(msg+leftchunkrecvd*CSIZE,CSIZE,MPI_CHAR,temp->right_child->process_id,leftchunkrecvd,MPI_COMM_WORLD,&sreq[count++]);
                    sent += 1;
                }
            }
            if (rightchunkrecvd != -1) {
                struct TreeNode* temp = rightTreeNode[rank];
                if(temp->left_child != NULL && temp->leftColor == turn) {
                    MPI_Isend(msg+rightchunkrecvd*CSIZE,CSIZE,MPI_CHAR,temp->left_child->process_id,rightchunkrecvd,MPI_COMM_WORLD,&sreq[count++]);
                    sent += 1;
                }
                else if (temp->right_child != NULL && temp->rightColor == turn) {
                    MPI_Isend(msg+rightchunkrecvd*CSIZE,CSIZE,MPI_CHAR,temp->right_child->process_id,rightchunkrecvd,MPI_COMM_WORLD,&sreq[count++]);
                    sent += 1;
                }
            }

            // left tree
            int willrecieve = 0;
            if ((leftTreeNode[rank]->parent->left_child->process_id == rank && leftTreeNode[rank]->parent->leftColor == turn)
                 || (leftTreeNode[rank]->parent->right_child->process_id == rank && leftTreeNode[rank]->parent->rightColor == turn)) {
                        j = leftchunkrecvd + 2;
                        int flag, MPI_Status st;
                        MPI_IPROBE(leftTreeNode[rank]->parent->process_id, j, MPI_COMM_WORLD, &flag, &st);
                        if (flag) {
                            MPI_Irecv(msg+j*CSIZE,CSIZE,MPI_CHAR,leftTreeNode[rank]->parent->process_id,j,MPI_COMM_WORLD,&req[j]);
                            willrecieve = 1;
                        }
            }
            else if ((rightTreeNode[rank]->parent->left_child->process_id == rank && rightTreeNode[rank]->parent->leftColor == turn)
                 || (rightTreeNode[rank]->parent->right_child->process_id == rank && rightTreeNode[rank]->parent->rightColor == turn)) {
                        j = rightchunkrecvd + 2;
                        int flag, MPI_Status st;
                        MPI_IPROBE(rightTreeNode[rank]->parent->process_id, j, MPI_COMM_WORLD, &flag, &st);
                        if (flag) {
                            MPI_Irecv(msg+j*CSIZE,CSIZE,MPI_CHAR,rightTreeNode[rank]->parent->process_id,j,MPI_COMM_WORLD,&req[j]);
                            willrecieve = 1
                        }
            }

            if (willrecieve == 1) {
                MPI_Waitany(CHUNK, req, &index, &stt);
	            	if(index == MPI_UNDEFINED) {
                		printf("Unexpected error!\n");
                		MPI_Abort(MPI_COMM_WORLD, 1);
            		}
                recieved += 1;
                if(index%2 == 0) {  // left tree
                    leftchunkrecvd = index;

                } else { //right tree
                    rightchunkrecvd = index;

                }
			}
            }
        }
    }
    MPI_Finalize();
    return 0;
}

