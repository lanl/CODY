//
//  QuadTree.h
//
/*
 *
 * When the update method is called the tree is refined and coarsened until each
 * cell meets the coarsening and refinement criteria which are determined at both
 * the global and local level through the inheritance from the base classes 
 * QuadTree and Application
 *
 * We developed a method for finding the face neighbors of a given node.  For
 * testing purposes we are currently outputting this data into a Gnuplot script
 * that can be used to show the neighbors.  Red nodes are north neighbors,
 * green nodes are south neighbors, magenta nodes are east neighbors,
 * and blue nodes are west neighbors
 *
 * Currently since we have implemented a pointer based quadtree, the getNeighbors
 * method requires a lot of pointer chasing.  There are more efficient indexing
 * schemes to implement the tree and then find the neighbors, but since I am
 * currently concerned with developing the overall interface not the
 * implementation, this is acceptable.
 */
//
//  Created by Alexandra Gendreau on 11/7/13.
//

#ifndef ____QuadTree__
#define ____QuadTree__

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <vector>
#include <cstdlib>
#include "Application.h"


/*
 * The QuadTree is composed of nodes
 */

struct Node{
    
public:
    
    Node * NEChild;//pointer to first quadrant child
    Node * NWChild;//pointer to second quadrant child
    Node * SWChild;//pointer to third quadrant child
    Node * SEChild;//pointer to fourth quadrant child
    Node * parent;//pointer to parent (used in neighbor finding and coarsening)
    Application * app; //what are we using this tree for
    double width; //how wide am I
    double height; //how tall am I
    double x; //where do I start x coord
    double y; //where do I start y coord
    int currentLevel;//level in tree where the node is
    int childType; //what child am I? 0 - NE, 1 - NW, 2-SW, 3-SE.  If root -1
    bool isLeaf;
    
    //constructor
    Node(double xStart, double yStart, double w, double h,
         Node * par, double level, int cType, Application * application){
     NEChild        = NULL;
     NWChild        = NULL;
     SWChild        = NULL;
     SEChild        = NULL;
     parent         = par;
     app            = application;
     x              = xStart;
     y              = yStart;
     width          = w;
     height         = h;
  
     currentLevel   = level;
     childType      = cType;
     isLeaf         = true;

     }
    
};



class QuadTree {
public:
    
    /*
     * Constructor that initializes the tree, takes in a pointer to root node,
     * the number of cells in the initial decomposition (must be a power of 4)
     * the maximum level
     */
    QuadTree(double x,double y, double width, double height,
             int numCells,int max,Application * app);
    
    /*
     * Destructor
     */
    ~QuadTree();
    
    /*
     * Refines and coarsens the quadtree until the desired refinement is reached
     */
    void update();
    
    /*
     * Refines and coarsens the node
     */
    void refineNode(Node * node);
    void coarsenNode(Node * node);
    
    /*
     * Returns false if the node cannot be coarsened or refined
     */
    virtual bool refine(Node * node);
    virtual bool coarsen(Node * node);

    
    /**********************DEBUGGING AND TREE INFO****************************/
    
    /* 
     * Returns pointers to all the leaf nodes
     */
    void findLeaves(std::vector<Node*>& leaves);
    
    /*
     * Method for getting the neighbors of a node in the tree
     */
    void getNeighbors(Node * node,
                      std::vector<std::vector<Node*> >& neighbors);
    
    /*
     * Method for getting dimensions of quad tree
     */
    std::vector<double> getDimensions();
    
    /*
     * Returns the total memory used to store the root node
     */
    int getSizeRoot();
    
    /* 
     * Getter and setter for max level of the tree
     */
    int getMaxLevel();
    void setMaxLevel(int level);
    
    /*
     * Getter and setter to turn update timing on/off
     */
    bool getTime();
    void setTime(bool t);
    
    /*
     * Returns the total amount of time spend in coarsening and refinement
     */
    double getTotalCoarsen();
    double getTotalRefine();
    
    /*
     * Counts the number of nodes in the tree
     */
    int countNodes();
    
    /*
     * Returns the total amount of memory used to store the quad tree
     */
    int storage();
    
    /*
     * Finds the leaf node that contains the particular x,y point
     */
    Node * findNode(float x,float y);
    
private:
    /*
     * Helper method for the constructor that constructs the initial spatial 
     * decomposition by inserting the correct number of nodes
     */
    void insert(Node * node, int levelsRemaining,int currentLevel);
    
    /*
     * Helper method for recursively destructing the tree
     */
    void destroyTree(Node * node);
    
    
    /*
     * Traverses the tree refining and coarsening the nodes
     */
    void checkCriteria(Node * node);
    
    /*
     * Refines the node as far as possible
     */
    void fullyRefine(Node * node);
    
/***************** DEBUGGING AND TREE INFO HELPERS ***************************/
    
    /*
     * Recursive helper methods for counting the nodes, printing the values,
     * and calculating the left Riemann summ, memory usage, finding the location
     * of an (x,y) pair, and the leaves of the tree.  These methods are private
     * in order to prevent the user from needing access to the root
     */
    
    int storageHelper(Node * node);
    int countNodesHelper(Node * node);
    Node * findNodeHelper(float x, float y, Node * node);
    void findLeavesHelper(std::vector<Node*>& leaves, Node * node);
    
    /*
     * Helper methods for finding the neighbors of a given node
     * There are four vectors with one vector
     * for the North Neighbors (first vector), South Neighbors (second vector),
     * East Neighbors (third vector), and West Neighbors (fourth vector)
     * These are face neighbors.  Each vector holds pointers to the neighbor nodes
     */

    std::vector<std::vector<Node*> > getParentsNeighbors(Node * parent);
    
    void getNeighborsSibs(Node * node, std::vector<Node*>& list,
                          int sibOneType, int sibTwoType);
    
    bool actualNeighbor(Node * node, int level, int nodeType);
    
    void processNeighbors(std::vector <Node*> pNeighbors,
                          std::vector<Node*>& list, int plevel, int nodeType);
    
/*********************** INSTANCE VARIABLES *********************************/
  
    Node * root;
    
    int maxLevel;
    bool time;
    double totalCoarsen;
    double totalRefine;
    
};

#endif /* defined(____QuadTree__) */