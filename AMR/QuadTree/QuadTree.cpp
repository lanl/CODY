//
//  QuadTree.cpp
//
//  See QuadTree.h for more detailed comments
//
//  Created by Alexandra Gendreau on 11/22/13.
//
//

#include "QuadTree.h"

using namespace std;

/*
 * Constructor
 */
QuadTree::QuadTree(double x, double y, double width, double height,
                   int numCells, int max, Application * app) {
    root            = new Node(x, y, width, height, NULL, 0, -1, app);
    int numLevels   = 0; //stores the number of levels in the tree
    int i           = numCells;
    while(i>>=1)
        ++numLevels;
 
    insert(root, (numLevels/2), 0);
    maxLevel        = max;
    time            = false;

}

/*
 * Constructor Helper
 */
void QuadTree::insert(Node * node, int levelsRemaining,int currentLevel){
    if(levelsRemaining == 1){
        //we reached the leaf level, insert nodes with non-zero value
        refineNode(node);
    }
    else if(levelsRemaining > 1){
        //we haven't reached the leaves yet so keep dividing
        int newLevel = node->currentLevel + 1;
        refineNode(node);
        
        insert(node->NEChild, levelsRemaining - 1, newLevel);
        insert(node->NWChild, levelsRemaining - 1, newLevel);
        insert(node->SWChild, levelsRemaining - 1, newLevel);
        insert(node->SEChild, levelsRemaining - 1, newLevel);
    }
}

/*
 * Destructor
 */
QuadTree::~QuadTree(){
    destroyTree(root);
}

/*
 * Recursive destructor helper
 */
void QuadTree::destroyTree(Node * node){
    if(node == NULL)
        return;
    else {
        destroyTree(node->NEChild);
        destroyTree(node->NWChild);
        destroyTree(node->SWChild);
        destroyTree(node->SEChild);
        
        delete node;
    }
}

/*
 * Updates the tree using refinement and coarsening criteria
 * Times the refinement and coarsening if time is true
 */
void QuadTree::update(){
    if(time){
        totalRefine     = 0;
        totalCoarsen    = 0;
    }
    checkCriteria(root);

}

/*
 * Recursive helper method for updating the tree
 */
void QuadTree::checkCriteria(Node * node){
    double start, finish;
    if(node->isLeaf){
        if(refine(node)){
            if(time)
                start       = clock();
            fullyRefine(node);
            if(time){
                finish      = clock();
                totalRefine = totalRefine + (double(finish-start)/CLOCKS_PER_SEC);
            }
        }
    }
    else if(coarsen(node)){
        if(time)
            start           = clock();
        coarsenNode(node);
        if(time){
            finish          = clock();
            totalCoarsen    = totalCoarsen + (double(finish-start)/CLOCKS_PER_SEC);
        }

    }
    else{
        checkCriteria(node->NEChild);
        checkCriteria(node->NWChild);
        checkCriteria(node->SWChild);
        checkCriteria(node->SEChild);
    }
}

/*
 * Returns a boolean determining if the node should be coarsened
 */
bool QuadTree::coarsen(Node * node){
        return node->app->coarsen(node->x,node->y,node->width,node->height);
}

/*
 * Returns a boolean determining if the node should be refined
 */
bool QuadTree::refine(Node * node){
    if(node->currentLevel == maxLevel){
              return false;
    }

    else {
        return node->app->refine(node->x,node->y,node->width,node->height);
    }
    
}

/*
 * Refines a leaf node by adding four children
 */
void QuadTree::refineNode(Node * node){
    double x        = node->x;
    double y        = node->y;
    double w        = (node->width)/2.0;
    double h        = (node->height)/2.0;
    int newLevel    = node->currentLevel +1;
    node->NEChild   = new Node((x+w),(y+h),w,h,node,newLevel,0,node->app);
    node->NWChild   = new Node(x,y+h,w,h,node,newLevel,1,node->app);
    node->SWChild   = new Node(x,y,w,h,node,newLevel,2,node->app);
    node->SEChild   = new Node(x+w,y,w,h,node,newLevel,3,node->app);
    node->isLeaf    = false;
}

/*
 * Coarsens a node by deleting all of its descendants
 */
void QuadTree::coarsenNode(Node * node){
    //cout<<"we tagged a node for coarsening"<<endl;
    destroyTree(node->SEChild);
    destroyTree(node->SWChild);
    destroyTree(node->NWChild);
    destroyTree(node->NEChild);
    
    //reset child pointers
    node->NEChild   = NULL;
    node->NWChild   = NULL;
    node->SWChild   = NULL;
    node->SEChild   = NULL;
    node->isLeaf    = true;
}

/*
 * Refines the nodes as far as necessary:
 *
 * To the maximum level or
 * The refinement criteria is no longer satisfied
 */
void QuadTree::fullyRefine(Node * node){
    
    refineNode(node);
    
    if(refine(node->NEChild))
        fullyRefine(node->NEChild);
    
    if(refine(node->NWChild))
        fullyRefine(node->NWChild);
    
    if(refine(node->SWChild))
        fullyRefine(node->SWChild);
    
    if(refine(node->SEChild))
        fullyRefine(node->SEChild);
    
}

/**********************DEBUGGING AND TREE INFO****************************/

int QuadTree::getSizeRoot(){
    return sizeof(Node);
}

int QuadTree::getMaxLevel(){
    return maxLevel;
}

void QuadTree::setMaxLevel(int level){
    maxLevel = level;
}

bool QuadTree::getTime(){
    return time;
}
void QuadTree::setTime(bool t){
    time = t;
}

double QuadTree::getTotalCoarsen(){
    return totalCoarsen;
}

double QuadTree::getTotalRefine(){
    return totalRefine;
}

//memory usage
int QuadTree::storage(){
    return storageHelper(root);
}

int QuadTree::storageHelper(Node * node){
    if(node->isLeaf)
        return sizeof(Node);
    else {
        return sizeof(Node)+
        storageHelper(node->NEChild)+
        storageHelper(node->NWChild)+
        storageHelper(node->SWChild)+
        storageHelper(node->SEChild);
    }
    
}

Node* QuadTree::findNode(float x, float y) {
    return findNodeHelper(x,y,root);
}

//find which node a given set of coordinates is located
Node* QuadTree::findNodeHelper(float x, float y, Node * node){
    if(node->isLeaf  && (x < (node->x + node->width)) &&
       (y < (node->y + node->height)))
        return node;
    else if ((x < node->x + (node->width)/2.0)){
        if(y < node->y + (node->height)/2.0)
            return findNodeHelper(x,y,node->SWChild);
        else
            return findNodeHelper(x,y,node->NWChild);
    }
    else {
        if(y < node->y + (node->height)/2.0)
            return findNodeHelper(x,y,node->SEChild);
        else
            return findNodeHelper(x,y,node->NEChild);
    }
    
}

void QuadTree::findLeaves(std::vector<Node*>& leaves){
    findLeavesHelper(leaves,root);
}

void QuadTree::findLeavesHelper(std::vector<Node*>& leaves,Node * node){
    if(node->isLeaf)
        leaves.push_back(node);
    else {
        findLeavesHelper(leaves,node->NEChild);
        findLeavesHelper(leaves,node->NWChild);
        findLeavesHelper(leaves,node->SWChild);
        findLeavesHelper(leaves,node->SEChild);
    }
}

//counts the number of nodes in tree recursively
int QuadTree::countNodes(){
    return countNodesHelper(root);
}

int QuadTree::countNodesHelper(Node * node){
    if(node->isLeaf)
        return 1;
    else{
        return 1+
        countNodesHelper(node->NEChild)+
        countNodesHelper(node->NWChild)+
        countNodesHelper(node->SWChild)+
        countNodesHelper(node->SEChild);
    }
}

std::vector<double> QuadTree::getDimensions(){
    vector<double> dims;
    dims.push_back(root->x);
    dims.push_back(root->y);
    dims.push_back(root->width);
    dims.push_back(root->height);
    return dims;
}


/************************* NEIGHBOR FINDING ******************************/

void QuadTree::getNeighbors(Node * node,
                            std::vector<std::vector<Node*> >& neighbors){
    if(node->parent == NULL){ //we're at the root which has no neighbors
        return;
    }
    
    vector<vector<Node*> > parentsNeighbors = getParentsNeighbors(node->parent);
    switch(node->childType){
        case 0://NEChild
            //cout<<"northeast child"<<endl;
            getNeighborsSibs(node->parent->NWChild,neighbors[3],0,3);
            getNeighborsSibs(node->parent->SEChild,neighbors[1],0,1);
            
            processNeighbors(parentsNeighbors[2],neighbors[2],
                             node->currentLevel-1,1);
            processNeighbors(parentsNeighbors[0],neighbors[0],
                             node->currentLevel-1,3);
            return;
            
        case 1: //NWChild
            //cout<<"northwest child"<<endl;
            getNeighborsSibs(node->parent->NEChild,neighbors[2],1,2);
            getNeighborsSibs(node->parent->SWChild,neighbors[1],0,1);
            
            processNeighbors(parentsNeighbors[3],neighbors[3],
                             node->currentLevel-1,0);
            processNeighbors(parentsNeighbors[0],neighbors[0],
                             node->currentLevel-1,2);
            return;
            
        case 2:  //SWChild
            //cout<<"southwest child"<<endl;
            getNeighborsSibs(node->parent->SEChild,neighbors[2],1,2);
            getNeighborsSibs(node->parent->NWChild,neighbors[0],2,3);
            
            processNeighbors(parentsNeighbors[3],neighbors[3],
                             node->currentLevel-1,3);
            processNeighbors(parentsNeighbors[1],neighbors[1],
                             node->currentLevel-1,1);
            return;
            
        case 3:  //SEChild
            //cout<<"southeast child"<<endl;
            getNeighborsSibs(node->parent->SWChild,neighbors[3],0,3);
            getNeighborsSibs(node->parent->NEChild,neighbors[0],2,3);
            
            processNeighbors(parentsNeighbors[2],neighbors[2],
                             node->currentLevel-1,2);
            processNeighbors(parentsNeighbors[1],neighbors[1],
                             node->currentLevel-1,0);
            return;
    }
    
    
}

/* Returns the neighbors of a node.
 * Used to calculate the parent nodes neighbors
 * allows recursive call without regenerating vector */
std::vector<std::vector<Node*> > QuadTree::getParentsNeighbors(Node * parent){
    //vector<vector<Node*> > parentsNeighbors;
    vector<Node*> vec;
    vector<vector<Node*> > parentsNeighbors (4, vec);
    getNeighbors(parent,parentsNeighbors);//get neighbors
    return parentsNeighbors;//return neighbors
}

void QuadTree::getNeighborsSibs(Node * node, std::vector<Node*>& list,
                                int sibOneType, int sibTwoType){
    if(node->isLeaf)
        list.push_back(node);
    else {
        if((sibOneType == 0) || (sibTwoType == 0)){
            getNeighborsSibs(node->NEChild,list,sibOneType,sibTwoType);
        }
        if(sibOneType == 1 || sibTwoType == 1){
            getNeighborsSibs(node->NWChild,list,sibOneType,sibTwoType);
        }
        if(sibOneType == 2 || sibTwoType == 2){
            getNeighborsSibs(node->SWChild,list,sibOneType,sibTwoType);
        }
        if(sibOneType == 3 || sibTwoType == 3) {
            getNeighborsSibs(node->SEChild,list,sibOneType,sibTwoType);
        }
        
    }
}

bool QuadTree::actualNeighbor(Node * node, int level, int nodeType){
    if(node->currentLevel == level)
        return node->childType == nodeType;
    else
        return actualNeighbor(node->parent,level,nodeType);
    
}

void QuadTree::processNeighbors(std::vector <Node*> pNeighbors,
                                std::vector<Node*>& list,int pLevel,
                                int nodeType){
    if(pNeighbors.size() == 0)
        return;
    
    else if(pNeighbors[0]->currentLevel <= pLevel)//less refined
        list.push_back(pNeighbors[0]);
    
    else{
        for(vector<Node*>::iterator vecIt = pNeighbors.begin();
            vecIt!=pNeighbors.end();vecIt++) {
            if((*vecIt)->currentLevel == pLevel+1 &&
               (*vecIt)->childType == nodeType)
                list.push_back((*vecIt)); //same level of refined
            else if(actualNeighbor((*vecIt),(pLevel+1),nodeType))//more refined have to check individually
                list.push_back((*vecIt));
        }
    }
}

/*int main(){
 
 }*/


