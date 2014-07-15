//
//  This is the main function that facilitates the drawing of the a QuadTree
//  object.  It initialize a pointer to a QuadTree object and uses openGL and
//  glut to draw the spatial decomposition of the quad tree and display
//  information such as neighbors using the treeRenderer class.
//  The treeRenderer class contains methods to draw the tree, neighbors, and
//  text using glut and the rectangular coordinates derived from the tree nodes
//
//  Created by Alexandra Gendreau on 1/7/14.
//
//

#include <cstdio>
#include "GLUT/glut.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <time.h>
//#include "QuadTree.h"
#include "OneLevel.h" //include the derived class
#include "Neighbor.h"
#include "treeRenderer.h"

//#include "Application.h"

using namespace std;

int NPIX; //resolution of image horizontal
int NPIY; //vertical
double leftX,leftY,width,height; //world info
int numCells = 16;//currently has to be a power of 4
bool displayNeighbors = false; //boolean for whether or not to display neighbors
bool coarsen = false; //boolean for whether or not to coarsen node
bool update = false; //determines if we are using the line node or graphics node
QuadTree* tree; //pointer to the QuadTree object we want to display
string text; //string for the text we display

Segment * seg; //pointer to line segment we are pushing through the space
double originalB;

//holds neighbor info, passed to treeRenderer
vector<vector<Rectangle> > neighborInfo;
//holds grid info, passed to treeRenderer
vector<Rectangle> gridInfo;

bool dynamic            = false;
double sleepTime        = 0;

bool updateDestructTest = true;
bool neighborTest       = true;
bool traversalTest      = true;
bool updateTest         = false;


/*
 * Converts (x,y) coordinate in pixel space to (x,y) coordinate in world spaces
 */
vector<double> treeLocal(double x, double y, int dimX,int dimY){
    vector<double> dims = tree->getDimensions();
    vector<double> act;
    double actX         = ((x*dims[2])/dimX)+dims[0];
    double actY         = (((dimY-y)*dims[3])/dimY)+dims[1];
    act.push_back(actX);
    act.push_back(actY);
    return act;
}

/* 
 * Converts a vector of Node objects to rectangular coordinates that can be
 * interpreted by the treeRenderer
 * Eventually a variation of this method might be in the quad tree class
 * instead of in the main function
 */
vector<Rectangle> convert(vector<Node*> nodes){
    vector<Rectangle> coords;
    for(vector<Node*>::iterator IT = nodes.begin();
        IT!=nodes.end();IT++){
        Rectangle rec((*IT)->x,(*IT)->y,(*IT)->width,(*IT)->height);
        coords.push_back(rec);
    }
    return coords;
}

/* 
 * Returns the rectangular coordinates of the neighbor of the node that was
 * clicked on in the window.
 * First we convert from pixel space to world space
 * Second we find the node that contains that world point
 * Third we find the neighbors of that node
 * Fourth we convert from nodes to rectangular coordinates to be drawn
 * Again this method or variation of it might be later placed in the
 * quad tree class
 */
vector<vector<Rectangle> > getNeighbors(double x,double y,
                                        int dimX,int dimY){
    vector<vector<Rectangle> > neighborCoords;
    vector<double> act = treeLocal(x,y,dimX,dimY);
    vector<Node*> vec;
    vector<vector<Node*> > neighbors (4, vec);
    //cout<<"neighbors:"<<act[0]<<","<<act[1]<<endl;
    Node * n = tree->findNode(act[0],act[1]);
    tree->getNeighbors(n,neighbors);
    for(vector<vector<Node*> >::iterator IT = neighbors.begin();
        IT!=neighbors.end();IT++){
        neighborCoords.push_back(convert(*IT));
    }
    return neighborCoords;
}

/* 
 * Returns the unfilled rectangles used to draw the grid
 * First finds all the leaves of the quad tree
 * Second convert leaves to rectangles to be drawn by treeRenderer
 * Again this method or variation of it might be later placed in the
 * quad tree class
 */
vector<Rectangle> getGrid(){
    vector<Node*> grid;
    tree->findLeaves(grid);
    return convert(grid);
}

/* 
 * Refines the node that was clicked in
 * First converts from pixel to world coordinates
 * Second finds the leaf node that contains the world coordinates
 * Third refine the node if possible (i.e. doesn't go over max levels)
 */
bool refineNode(double x, double y, int dimX,int dimY){
    vector<double> act = treeLocal(x,y,dimX, dimY);
    Node * node = tree->findNode(act[0],act[1]);
    if(tree->refine(node)){
        tree->refineNode(node);
        return true;
    }
    else
        return false;
}

/* 
 * Coarsens the node that was clicked in up one level
 * First converts from pixel to world coordinates
 * Second finds the leaf node that contains the world coordinates
 * Third coarsen the node if possible (i.e. isn't the root or one or more of
 * its siblings aren't leaves)
 */
bool coarsenNode(double x,double y,int dimX,int dimY){
    vector<double> act  = treeLocal(x,y,dimX,dimY);
    Node * node         = tree->findNode(act[0],act[1]);
    if(tree->coarsen(node)){
        tree->coarsenNode(node->parent);
        return true;
    }
    else
        return false;
}

/*
 * Updates the spatial decomposition based on the current position of the line
 */
void updateTree(){
    tree->update();
}

/*
 * The function that calls that calls the draw functions in the treeRenderer class
 */
void display(){
    /*  clear all pixels  */
    
    glClear (GL_COLOR_BUFFER_BIT);
    glColor3f (1.0, 1.0, 1.0);
    if(dynamic)
        updateTree();
    gridInfo = getGrid();
    
    
    treeRenderer::drawString(GLUT_BITMAP_HELVETICA_18, text, leftX, leftY, 0);
    
    if(displayNeighbors)
        treeRenderer::displayHelperNeighbors(neighborInfo);
    
    treeRenderer::displayHelperTree(gridInfo);
    
    if(update || dynamic)
        treeRenderer::drawLine(seg->getx0(),seg->gety0(),
                               seg->getx1(),seg->gety1()); //draws line
    
    glFlush ();
    
    if(dynamic){
        if(tree->countNodes()==1) //the line is no longer in our space, reset
            seg->reset(originalB);
        else
            seg->translate(1.0/NPIX,1.0/NPIY);
    }
    
    
    glutSwapBuffers();
}

/*
 * Initialize the window
 */
void init (){
    /*  select clearing (background) color       */
    glClearColor (0.0, 0.0, 0.0, 0.0);
    
    /*  initialize viewing values  */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(leftX, leftX+width, leftY, leftY+height);//sets window size
    
}

/*
 * Facilitates interactions between the keyboard and viewing window
 */
void keyboard(unsigned char key, int x, int y){
    switch (key) {
        case 'c':{//mouse clicks cause coarsening
            if(!coarsen)
                coarsen = true;
            else
                coarsen = false;
            break;
        }
        case 'n': {//mouse click displays neighbors
            if(!displayNeighbors)
                displayNeighbors = true;
            else
                displayNeighbors = false;
            break;
        }
        case 'm':{ //use the complete update method instead of clicking
            if(update)
                update = false;
            else
                update = true;
            break;
        }
        case 't':{ //translate the segment so it progresses through the space
            if(tree->countNodes()==1) //the line is no longer in our space, reset
                seg->reset(originalB);
            else
                seg->translate(0.25,0.25);
            
            glutPostRedisplay();
            break;
        }
        case 'd':{
            if(dynamic){
                dynamic     = false;
                sleepTime   = 0;
            }
            else {
                dynamic     = true;
                sleepTime   = 0.2;
            }
            break;
        }
        case 27: //escape key
            delete tree; // free up the memory, is the correct?
            exit(0); //exit
            break;
        default:
            break;
    }
}

/*
 * Facilitates interactions between the mouse and the viewing window
 */
void mouse(int button, int state, int x, int y){
    // Respond to mouse button presses.
    if (state == GLUT_DOWN){
        if(update){
            updateTree();
            //cout<<"we finished updating"<<endl;
            gridInfo        = getGrid();
            
        }
        if(displayNeighbors){//if the 'n' key was pressed, display neighbors
            neighborInfo    = getNeighbors(x,y,NPIX,NPIY);
            text            = "";
        }
        if(!update && !displayNeighbors){//refine and coarsen using mouse clicks
            if(coarsen){//if the 'c' key was pressed, coarsen the node
                if(!coarsenNode(x,y,NPIX,NPIY))
                    text    = "Cannot Coarsen: at root or siblings not leaves";
                else {
                    gridInfo = getGrid();
                    text    = "";
                }
            }
            else {//if neither key was pressed and we are not using the complete
                //update, refine the node
                if(!refineNode(x,y,NPIX,NPIY))
                    text        = "Cannot Refine: reached max level";
                else {
                    gridInfo    = getGrid();
                    text        = "";
                }
            }
        }
        glutPostRedisplay();//display changes
    }
}

/*
 * Reshapes the window if the user adjusts the size
 */
void reshape(int w, int h){
    NPIX = w;
    NPIY = h;
    glViewport(0, 0, (GLsizei) NPIX, (GLsizei) NPIY);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(leftX, leftX+width, leftY, leftY+height);
    
}


/* 
 * OpenGL idle function that sets the frame rate
 * Makes the flow visualization easier to see
 */
void idle(){
    
    sleep(sleepTime);
    glutPostRedisplay();
    
}

/* main
 * Initializes:
 *  the QuadTree object
 *  the window size
 *  glut
 * Runs glut main loop
 */
int main(int argc, char* argv[]){
    //QuadTree and window size initialization
    if(string(argv[1]) == "graphics"){ //use graphcis nodes
        cout<<"here"<<endl;
        if(string(argv[2]) == "step") {
            Interaction * app   = new Interaction();
            update              = false;
            dynamic             = false;
            tree                = new OneLevel(-4.0,-4.0,4.0,4.0,16,4,app);
        }
        else if(string(argv[2]) == "line") { //use line nodes
            update              = true;
            originalB           = -7.0;
            seg                 = new Segment(-4.0,0.0,-1.0,originalB);
            Line * app2         = new Line(seg);
            tree                = new QuadTree(-4.0,-4.0,4.0,4.0,16,4,app2);
        }
        else if(string(argv[2]) == "neighbor") { //neighbor restriction
            Interaction * app = new Interaction();
            update              = false;
            dynamic             = false;
            tree                = new Neighbor(-4.0,-4.0,4.0,4.0,16,4,app);
        }
        
        vector<double> dims     = tree->getDimensions();
        leftX                   = dims[0];
        leftY                   = dims[1];
        width                   = dims[2];
        height                  = dims[3];
        NPIX = NPIY             = 512;
        gridInfo                = getGrid();
        
        //glut initialization
        glutInit(&argc, argv);
        glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
        glutInitWindowSize (NPIX,NPIX);
        glutInitWindowPosition (100, 100);
        glutCreateWindow ("QuadTree");
        glutDisplayFunc(display);
        glutKeyboardFunc(keyboard);
        glutMouseFunc(mouse);
        glutIdleFunc(idle);
        glutReshapeFunc(reshape);
        
        init ();
        glutMainLoop();
    }
    else if(string(argv[1]) == "run"){ //we can do some of our testing here
    //outputs csv files of testing results
        originalB               = -7.0;
        seg                     = new Segment(-4.0,0.0,-1.0,originalB);
        seg->translate(1.5,1.5);
        Line * app3             = new Line(seg);
        clock_t start,finish;
        
        if(updateDestructTest) {
            //Tests the amount of time for the constructor and destructor as
            //as well as memory for different sizes trees
            ofstream file;
            file.open("updateDestruct.csv");
            file <<"nodes,memory,update,destruct \n";
            
            for(int i=0;i<10;i++){//construct trees with a line through the center of different maximum depths, also neighbors
                double averageTimeUpdate = 0.0;
                double averageTimeDestructor = 0.0;
                int nodes;
                int memory;
                for(int j=0;j<10;j++){
                    tree    = new QuadTree(-4.0,-4.0,4.0,4.0,1,i,app3);
                    start   = clock();
                    updateTree();
                    finish  = clock();
                    nodes   = tree->countNodes();
                    memory  = tree->storage();
                    averageTimeUpdate += (double(finish-start)/CLOCKS_PER_SEC);
                
                    start   = clock();
                    delete tree;
                    finish  = clock();
                    averageTimeDestructor += (double(finish-start)/CLOCKS_PER_SEC);
                }
                averageTimeUpdate       = averageTimeUpdate/10.0;
                averageTimeDestructor   = averageTimeDestructor/10.0;
                if(file.is_open())
                    file <<nodes<<","<<memory<<","<<averageTimeUpdate<<","<<
                    averageTimeDestructor<<"\n";
                else
                    cout<<"FILE ERROR"<<endl;
            }
            file.close();
            
        }
        
        if(neighborTest){
            // Tests the amount of time to find all the neighbors of all the
            // nodes for trees of different sizes
            tree = new QuadTree(-4.0,-4.0,4.0,4.0,1,0,app3);
            ofstream file;
            file.open("neighborTest.csv");
            file <<"leaves,time,neighbors \n";
            for(int i=0;i<10;i++){//construct trees with a line through the center of different maximum depths, also neighbors
                double averageTimeNeighbors = 0.0;
                int nodes;
                int numberOfNeighbors = 0;
                for(int j=0;j<10;j++){
                    vector<Node *> leaves;
                    tree    = new QuadTree(-4.0,-4.0,4.0,4.0,1,i,app3);
                    updateTree();
                    nodes   = tree->countNodes();
                    tree->findLeaves(leaves);
                    nodes   = leaves.size();
                    start   = clock();
                    for(int m = 0;m<leaves.size();m++) {
                        vector<Node*> vec;
                        vector<vector<Node*> > neighbors (4, vec);
                        tree->getNeighbors(leaves[m],neighbors);
                        numberOfNeighbors += neighbors[0].size();
                        numberOfNeighbors += neighbors[1].size();
                        numberOfNeighbors += neighbors[2].size();
                        numberOfNeighbors += neighbors[3].size();
                    }
                    finish = clock();
                    averageTimeNeighbors += (double(finish-start)/CLOCKS_PER_SEC);
                    delete tree;
                }
                averageTimeNeighbors    = averageTimeNeighbors/10.0;
                numberOfNeighbors       = numberOfNeighbors/10;
                if(file.is_open())
                    file <<nodes<<","<<averageTimeNeighbors<<","<<numberOfNeighbors<<"\n";
                else
                    cout<<"FILE ERROR"<<endl;
            }
            file.close();
        }
        
        if(traversalTest){
            // Tests the total traversal time with trees of varying sizes
            tree = new QuadTree(-4.0,-4.0,4.0,4.0,1,0,app3);
            ofstream file;
            file.open("traversalTest.csv");
            file <<"nodes,traversal \n";
            for(int i=0;i<10;i++){//construct trees with a line through the center of different maximum depths, also neighbors
                double averageTimeTraversal = 0.0;
                int nodes;
                for(int j=0;j<10;j++){
                    tree    = new QuadTree(-4.0,-4.0,4.0,4.0,1,i,app3);
                    updateTree();
                    start   = clock();
                    nodes   = tree->countNodes();
                    finish  = clock();
                    averageTimeTraversal += (double(finish-start)/CLOCKS_PER_SEC);
                    delete tree;
                }
                averageTimeTraversal = averageTimeTraversal/10.0;
                if(file.is_open())
                    file <<nodes<<","<<averageTimeTraversal<<"\n";
                else
                    cout<<"FILE ERROR"<<endl;
            }
            file.close();
        }
        
        if(updateTest){
            
            //testing the percentage of time spend on coarsening and refinement
            ofstream file;
            file.open("updateTest.csv");
            file <<"initial nodes,final nodes,refine,coarsen,update \n";
            for(int k = 10;k<23;k++) {
                seg->translate(-2.0,-2.0);
                tree = new QuadTree(-4.0,-4.0,4.0,4.0,1,k,app3);
                tree->setTime(true);
                double refineTime   = 0.0;
                double coarsenTime  = 0.0;
                double updateTime   = 0.0;
                int initialNodes;
                int finalNodes;
                initialNodes        = tree->countNodes();
                start = clock();
                updateTree();
                finish = clock();
                updateTime          = (double(finish-start)/CLOCKS_PER_SEC);
                finalNodes          = tree->countNodes();
                refineTime          = tree->getTotalRefine();
                coarsenTime         = tree->getTotalCoarsen();
                seg->translate(0.25,0.25);
                initialNodes        = tree->countNodes();
                start               = clock();
                updateTree();
                finish              = clock();
                updateTime          = (double(finish-start)/CLOCKS_PER_SEC);
                finalNodes          = tree->countNodes();
                refineTime          = tree->getTotalRefine();
                coarsenTime         = tree->getTotalCoarsen();
                seg->translate(0.25,0.25);
                if(file.is_open())
                    file <<initialNodes<<","<<finalNodes<<","<<refineTime<<","<<
                    coarsenTime<<","<<updateTime<<"\n";
                else
                    cout<<"FILE ERROR"<<endl;
                while(tree->countNodes() > 1){
                    initialNodes    = tree->countNodes();
                    start           = clock();
                    updateTree();
                    finish          = clock();
                    updateTime      = (double(finish-start)/CLOCKS_PER_SEC);
                    finalNodes      = tree->countNodes();
                    refineTime      = tree->getTotalRefine();
                    coarsenTime     = tree->getTotalCoarsen();
                    seg->translate(0.25,0.25);
                    if(file.is_open())
                        file <<initialNodes<<","<<finalNodes<<","<<refineTime<<","<<
                        coarsenTime<<","<<updateTime<<"\n";
                    else
                        cout<<"FILE ERROR"<<endl;
                }
                file<<"\n";
                delete tree;
                seg->translate(-2.25,-2.25);
            }
            file.close();
       
        }
    }
    cout<<"we exited normally"<<endl;
    return 0;
}
