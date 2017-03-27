#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
	
int e; ///A global variable to hold the number of nodes in the tree which is a user-input dynamic value.
int *a; /// This global array holds the x coordinates of the nodes in the tree. 
int *b;/// This global array holds the y coordinates of the nodes in the tree.
int *stat; /// This global array hold the status of the nodes of the tree.
int gh; /// This variable is given the value maxh+1: This is used to invert the y coordinate which is assigned in an order that places root at the bottom and leaves at the top in the TR algorithm.
int c=0; /// A global variable to hold the index to the arrays a[], b[] and stat[], which contain the x, y and z coordinates of each node traversed in pre-order.
int max_x=0, max_y=0; /** Variables to hold the maximum of x and y values calculated from TR algorithm for a given tree. This is used for scaling the tree appropriately to the size of the monitor.*/
int min_x=1366,min_y=1366; /** Variables to hold the minimum of x and y values calculated from TR algorithm for a given tree. This is used for scaling the tree appropriately to the size of the monitor.*/
int final[3][1]; /** This column vector is used to calculate the homogenised transformation of x, y coordinates held in pt[3][1] to obtain the scaled and translated points. */
int pt[3][1]; /** This column vector holds the homogenised x coordinate of a point calculated through the TR algorithm, and the homogenised y coordinate after inversion (root at the top, leaves at the bottom)*/
int root[3][1]; /** This column vector hold the homogenised x and y coordinates of the root of the tree.*/
int temp[3][1]; /**This column vector hold the temporary values while calculating the transformation of each homogenised point corresponding to the nodes of the tree.*/
int Trans1[3][3] = {0};/** This is the translation matrix used to translate the given coordinates to the origin, initialized to 0.*/
int Scale[3][3] = {0}; /** This is the scaling matrix to be applied to the coordinates at the origin, initialized to 0.*/
int Trans2[3][3]={{1, 0,50}, {0, 1, 50}, {0, 0, 1}}; /**This is the translation matrix applied to points scaled at the origin, to place the tree in the viewport.*/
clock_t start, end; /// Clock variables to hold the number of clock cycles since the start and at the end of the program, used to calculate the runtime of the program.
double cpu_time; /// Variable used to hold the total runtime of the program.

float sx; /**
sy holds the scaling ratio along y axis.
*/
float sy;/**
*sx holds the scaling ratio along x axis.
*/


typedef struct Node{
	struct Node* parent; ///The pointer to the parent of the current node.
	struct Node *left; /// The pointer to the left child of the current node.
	struct Node *right; /// The pointer to the right child of the current node.
	int x; /// Variable to hold the x value of the Node.
	int y; /// Variable to hold the y coordinate of the Node.
	int height; /// This variable holds the height of the current node as it is used in implementing the TR algorithm. This is closely related to the variable y. 
	int status; /// Boolean value to indicate if the current node is a leaf (1) or a non-leaf node (0). 
	int modifier; /// Value of modifier used to adjust the positions of the nodes.
	int value;	///Holds the value of the node, used to determine the position of the node in the tree.
	}node;
node *God; ///Global pointer to the root of the tree.

void init(void){/**
This function initializes the graphics window. It clears the background colour to white, sets the matrix mode to projection and enables 2-D drawing by projecting and clipping the viewport to the screen size of 1366x768 pixels.
*/
  glClearColor(1,1,1,0);
  glMatrixMode(GL_PROJECTION);
  gluOrtho2D(0.0,1366.0,0.0,768.0); 
}

void setPixel(GLint x,GLint y){/**
This function sets the drawing mode to points, so as to enable pixel by pixel colouring, and sets the pixel specified by the arguments as x and y coordinates.
*/
  glBegin(GL_POINTS);
  glVertex2i(x,y);
  glEnd();
  glFlush();
}

void Circle(int x1,int y1){/**
This function uses the midpoint circle drawing algorithm to draw circles pixel by pixel, at node coordinates determined by TR algorithm. 
*/
	int xCenter=x1; /// The x coordinate of the center of the circle to be drawn. 
	int yCenter=y1; /// The y coordinate of the center of the circle to be drawn.
	int r=5; /// The radius of the circle to be drawn. We have taken a constant radius of 5 for drawing trees in this program.
	int x=0; 
	int y=r;
	int p = 3/2 - r;
	while(x<=y){
		setPixel(xCenter+x,yCenter+y);
		setPixel(xCenter+y,yCenter+x);
    	setPixel(xCenter-x,yCenter+y);
    	setPixel(xCenter+y,yCenter-x);
    	setPixel(xCenter-x,yCenter-y);
    	setPixel(xCenter-y,yCenter-x);
    	setPixel(xCenter+x,yCenter-y);
    	setPixel(xCenter-y,yCenter+x);
		if (p<0)
			p += (2*x)+3;
		else{
 			p += (2*(x-y))+5;
 			y -= 1;
    	}
    	x++;
  	}
}


void line (node *nodec,int x1,int y1){/**
This function takes in the pointer to a node, along with transformed x and y coordinate of the node's parent as its arguments.The coordinates of the node undergo translation and scaling and translation again, to determine the transformed coordinates of the node. A circle is then drawn using Circle() at the specified node's coordinates, after which Bresenham's line drawing algorithm is used to connect the transformed coordinates of the node to the coordinates of the parent node passed as arguments to this function. This is done in a pre-order fashion.
*/
	if(nodec == NULL)
		return;

	int x0=nodec->x; /// The x coordinate of the node before transformation.
	int y0=gh-nodec->y; /// The y coordinate of the node before transformation/
	int i, j, k; /// Indices
	
	
	pt[0][0] = x0;
	pt[1][0] = y0;
	pt[2][0] = 1;
	
	for(i=0;i<3;i++)
		final[i][0] = 0;

 	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			final[i][0]+=Trans1[i][j]*pt[j][0];			///Translating to origin
		}
	}
	pt[0][0] = final[0][0];
	pt[1][0] = final[1][0];
	pt[2][0] = final[2][0];

	for(i=0;i<3;i++)
		final[i][0] = 0;


	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			final[i][0]+=Scale[i][j]*pt[j][0];			/// Scaling by a factor of sx and sy along x and y respectively, at the origin.
		}
	}
	pt[0][0] = final[0][0];
	pt[1][0] = final[1][0];
	pt[2][0] = final[2][0];


	for(i=0;i<3;i++)
		final[i][0] = 0;
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			final[i][0]+=Trans2[i][j]*pt[j][0];			///Translating to (50, 50)
		}
	}
	x0 = final[0][0];
	y0 = final[1][0];

	int currentx = x0; /// The transformed x coordinate of the node.
	int currenty = y0; /// The transformed y coordinate of the node.

	if(nodec->status == 1){
		glColor3f(0.134, 0.545,0.134); //Green leaf
	}
	Circle(x0, y0);
	glColor3f(1, 0, 0.5); //fuchsia
	int tempx, tempy; /// Temporary variables for swapping
	if(x1<x0){
		tempx = x1;
		x1=x0;
		x0=tempx;
		tempy = y1;
		y1=y0;
		y0 = tempy;
		}

	int dx=x1-x0; /// x increment used in Bresenham's line drawing algorithm.
	int dy=y1-y0; /// y increment used in Bresenham's line drawing algorithm.
	float m=(float)dy/dx;
	if(m>0 && m<=1){
   		int d=2*dy-dx; /// The decision variable used to choose between two adjacent pixels while drawing the line.
   		int ie=2*dy; /// The increment for "eastern" pixel.
   		int ine=2*(dy-dx); /// The increment for "north-eastern" pixel.
   		setPixel(x0,y0);
   		while(x0<x1){
   			if(d<=0){
      				d=d+ie;
      			}
   			else{
      				d=d+ine;
      				y0++;
     			}
   			x0++;
   			setPixel(x0,y0);
   		}
	}	

	else if(m>1){
		int d=(-2*dx+dy);
		int ie=(-2*dx);
		int ine=2*(dy-dx);
		setPixel(x0,y0);
		while(x0<x1){
   			if(d>0){
    				d=d+ie;
    			}
   			else{
				d=d+ine;
				x0++;
     			}
   			y0++;
			setPixel(x0,y0); 
		}

	} 

	else if((m<0)&&(m>=-1)){
		int d=-(2*dy+dx);
		int ie=-2*dy;
		int ine=-2*(dy+dx);
		setPixel(x0,y0);
		while(x0<x1){
			if(d<=0){
				d=d+ie;
			}
			else{
				d=d+ine;
				y0--;
     			}
   			x0++;
   			setPixel(x0,y0); 
     		} 
   	}


  	else{
  		int d=-(2*dx+dy);
  		int ie=-2*dx;
   		int ine=-2*(dy+dx);
   		setPixel(x0,y0);
   		while(x0<x1){
      			if(d>0){
         			d=d+ie;
         		}
      			else{
         			d=d+ine;
         			x0++;
         		}
      		y0--;
      		setPixel(x0,y0); 
     		}  
  	}  
/**Recursive calls to the children nodes.*/
	line(nodec->left,currentx,currenty);
	line(nodec->right,currentx, currenty);
	return;
  	}


void drawRoot(int x1,int y1){/**
This function determines the transformed coordinates of the root of the tree and calls line(). This is done because line() expects tranformed coordinates of the parent node and the current node as arguments. This also draws a circle at the tranformed coordinates of the root so as to initialize the recursion.
*/

///Initializing the column vectors holding coordinates of the node.
	root[0][0] = 0;
	root[1][0] = 0;
	root[2][0] = 0;
	temp[0][0] = x1;
	temp[1][0] = y1;
	temp[2][0] = 1;
	int i,j;
	
 	for(i = 0;i<3;i++){
		for(j=0;j<3;j++){
			root[i][0] += Trans1[i][j] * temp[j][0];		/// Translating to origin.
		}
	}
	
	temp[0][0] = root[0][0];
	temp[1][0] = root[1][0];	
	temp[2][0] = root[2][0];
	
	root[0][0] = 0;
	root[1][0] = 0;
	root[2][0] = 0;

	for(i = 0;i<3;i++){
		for(j=0;j<3;j++){
			root[i][0] += Scale[i][j] * temp[j][0];		/// Scaling by a factor of sx and sy.
		}
	}
	temp[0][0] = root[0][0];
	temp[1][0] = root[1][0];	
	temp[2][0] = root[2][0];

	root[0][0] = 0;
	root[1][0] = 0;
	root[2][0] = 0;

	for(i = 0;i<3;i++){
		for(j=0;j<3;j++){
			root[i][0] += Trans2[i][j] * temp[j][0];		/// Translating the root to (50, 50)
		}
	}

/**
Implementing the midpoint circle drawing algorithm to draw a circle at the coordinates of the root.
*/			
	
	int xCenter = root[0][0];
	int yCenter = root[1][0];
	int r=5;
	int x=0,y=r;
	int p = 3/2 - r;
	while(x<=y){
		setPixel(xCenter+x,yCenter+y);
		setPixel(xCenter+y,yCenter+x);
    		setPixel(xCenter-x,yCenter+y);
    		setPixel(xCenter+y,yCenter-x);
    		setPixel(xCenter-x,yCenter-y);
    		setPixel(xCenter-y,yCenter-x);
    		setPixel(xCenter+x,yCenter-y);
    		setPixel(xCenter-y,yCenter+x);
		if (p<0)
			p += (2*x)+3;
		else{
 			p += (2*(x-y))+5;
 			y -= 1;
    		}
    		x++;
  	}
	glColor3f(1, 0, 0.5);
	line(God->left,xCenter,yCenter);
	line(God->right,xCenter, yCenter);

}
void display(){	/**
This function sets the size of the point to be drawn, calculates the minimum and maximum of the untransformed coordinates by traversing through all the nodes to proportionately place the tree on the viewport. The scaling and translation matrices are thus calculated from the minimummand maximum coordinates. This function also calls drawRoot() from which the drawing of the tree begins.
*/
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(3.2); 				///Setting the point size to 3.2
	glColor3f(1, 0, 0.5);			/// Fuchsia
	for(c=0;c<e;c++){
		if(a[c]<min_x)  /// for calculating max and min of the x and y coordinates.
			min_x=a[c];
		
		if(a[c]>max_x){
			max_x=a[c];     	
		}
		if(b[c]<min_y){
			min_y=b[c];
		}
		if(b[c]>max_y){
			max_y=b[c];     	
		}
    	}
	sx=1316/(max_x-min_x+2);		///Scaling factor along x
 	sy=700/(max_y-min_y+2);			///Scaling factor along y

	Scale[0][0] = sx;
	Scale[1][1] = sy;
	Scale[2][2] = 1;
	Trans1[0][0] = 1;
	Trans1[1][1] = 1;
	Trans1[0][2] = -min_x;
	Trans1[1][2] = -min_y;
	Trans1[2][2] = 1;
	glColor3f(0.435, 0.306, 0.216); //Brown root
	drawRoot(God->x, gh-(God->y));
	glFlush();
	end = clock(); /** 
After the drawing is completed, the code reaches this point where we call clock() to calculate the total time taken for drawing the tree.
*/
	cpu_time = (double)(end - start)/CLOCKS_PER_SEC;
	printf("Time: %lf\n", cpu_time);
}




node *createNode(int v){/**
This function is used to create the nodes for the tree dynamically. It assigns the node with a value, and a position with respect to the other nodes of the tree. This also initializes a number of parameters with a default value that could be changed during the run of the program to determine the exact coordinates of the node in the tree. 
*/
	node *node1 = (node*)malloc(sizeof(node));
	node1->value = v;
	node1->left =NULL;
	node1->right = NULL;
	node1->status = 1; //leaf when inserted.
	node1->height = 1;
	node1->parent = NULL;
	node1->x = 0;
	node1->y = 0; 
	return node1;
}

void postOrder(node *node1, int nextpos[], int modifier[]){/**
The function uses a post order traversal to determine provisional x coordinates based on the next available position at each level maintained in the next_pos array. The modifier is calculated to maintain ample distance between the left and right children at each node. At each level, it is calculated and stored in modifier array.
*/
	if(node1 == NULL)
		return;
	postOrder(node1->left, nextpos, modifier);
	postOrder(node1->right, nextpos, modifier);
	int place;
	if(node1->status == 1)
		place = nextpos[(node1->height)-1];
	else if(node1->left == NULL)
		place = ((node1->right)->x - 1);
	else if(node1->right == NULL)
		place = ((node1->left)->x + 1);
	else place = ((node1->left)->x +(node1->right)->x)/2;

	if(modifier[(node1->height)-1] <(nextpos[(node1->height)-1] - place))
		modifier[(node1->height)-1] = nextpos[(node1->height)-1] - place;

	if(node1->status == 1)
		node1->x =place;
	else
		node1->x = place+modifier[(node1->height)-1];
	nextpos[(node1->height)-1] = (node1->x) + 2;
	node1->modifier = modifier[(node1->height)-1];

	return;
}

void preOrder(node *node1){/**
The modifiers are added cumulatively from top to bottom of the tree using a pre-oreder traversal and the final coordinates of the nodes are determined
*/
	if(node1 == NULL)
		return;
	node1->y = 2*(node1->height) - 1;
	if(node1->parent != NULL){
		node1->x = (node1->x) + ((node1->parent)->modifier);
		node1->modifier = node1->modifier + (node1->parent)->modifier;
	}
	preOrder(node1->left);
	preOrder(node1->right);

	return;
}

void dcir(node *node1){ /** 
*This function is used to store the x and y coordinates, and status (leaf/non-leaf) of the nodes in global arrays. The function traverses through the nodes in pre-order. 
*/
	if(node1==NULL)
		return;

	a[c]=node1->x; /// This array stores the x corrdinates calculated by implementing the TR algorithm.
	b[c]=gh-(node1->y); /// This array stores the y coordinates calculated by implementing the TR algorithm.
	stat[c] = (node1->status);	/// This array stores the leaf/non-leaf status of each node.
    c++;
	printf("%d-> %d %d\n", node1->value, node1->x,node1->y);
    dcir(node1->left);
    dcir(node1->right);
    return;
    }


int main(int argc,char **argv){

	start = clock();
	int h; //Height of the tree for each node inserted: Calculated for each node to keep track of the max depth of the tree
	int val; // Values to scan from the user to build a tree
	int num; // Number of elements in the tree
	int i; // index
	printf("Enter the number of nodes in the tree:\n");
	scanf("%d", &num); //User inputs the number of elements	
	if(num == 0) return 0;
	printf("Would you like to enter %d inputs yourself (press 0) or would you like %d randomly generated inputs? (press 1)\n", num, num);
	int option; /** This variable takes an option from the user: 0 for user input nodes and 1 for randomly generated values for the nodes*/
	scanf("%d", &option);
	node *temp, *current; // Two nodes to keep a track of nodesto make comparisons of values  while building the tree
	time_t t;	
	srand((unsigned)time(&t));
	if(option == 1){
		val = rand()%1000;
		printf("Node 1: %d\n",val);	
	}
	else{
		printf("Node 1: ");
		scanf("%d", &val);
	}	
	node *root = createNode(val);
	root->parent=NULL;
	current = root;
	e=num;
	int maxh=1; //This variable contains the maximum depth of the tree required to draw a tidy tree
	/*Scanning inputs from the user and building the binary tree*/
	for(i=1;i<num;i++){
		h=1;		//Each node potentially begins at height 1 and as its position is determined, h is incremented.
		if(option == 1){
			val = rand()%1000;
			printf("Node %d: %d\n", i+1, val);		
		}
		else scanf("%d", &val);
		temp = createNode(val);
		while(current != temp){
			if(current->value < temp->value){
				if(current->right != NULL){
					current=current->right;
					h++;
				}
				else{
                    current->status =0; ///As soon as we input a new node to the right or left of a node, it becomes parent to the new node, hence it is not a leaf anymore.
					current->right = temp;
					temp->parent = current;
					current = temp;
					h++;
				}
			}
			else{
				if(current->left !=NULL){
					current = current->left;
					h++;
				    }
				else{	
					current->status = 0; //becomes parent to the new node, hence leaf no more.
					current->left = temp;
					temp->parent = current;
					current = temp;
					h++;
				   }
		      }
		  }

		current->height = h;
		if(maxh<h)
			maxh = h;
		current = root;
	}

	int nextpos[maxh]; /** 
This array holds the next available position at each level in the tree. This gets updated dynamically as we place the nodes in a post order traversal. 
*/
	int modifier[maxh]; /**
This array incrementally adds the x separation from the parent, so as to appropriately place the left child to the left and the right child to the right of a node.
*/

	for(i=0;i<maxh;i++){
		nextpos[i] = 1;
		modifier[i] = 0;
	}
	postOrder(root, nextpos, modifier);
	preOrder(root);

	God=root; //for the line algo
	gh=2*maxh;
	printf("\n");
 	a = (int *)malloc(sizeof(int)*num); //set the global variables
	b = (int *)malloc(sizeof(int)*num);
	stat = (int *)malloc(sizeof(int)*num);
	dcir(root); //pass the coordinates of the nodes to the arrays

 //Initialising OpenGL
   	glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(1366,768);
    glViewport(0,0,1366,768);
    glutCreateWindow("Tidy Drawings- Wetherell and Shannon");
    init();
    glutDisplayFunc(display);
    glutMainLoop();
	return 0;
}
