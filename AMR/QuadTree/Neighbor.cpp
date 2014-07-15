//
//  Neighbor.cpp
//  
//
//  Created by Alexandra Gendreau on 2/14/14.
//
//

#include "Neighbor.h"

using namespace std;

Neighbor::Neighbor(double x,double y, double width, double height,
                   int numCells, int max, Application * app) : QuadTree(x,y,width,height, numCells, max, app){
}

Neighbor::~Neighbor(){}

bool Neighbor::refine(Node * node){

    return QuadTree::refine(node) && levelPassed(node);
}

bool Neighbor::coarsen(Node * node){
    return QuadTree::coarsen(node) && node->parent->NEChild->isLeaf &&
    node->parent->NWChild->isLeaf &&
    node->parent->SWChild->isLeaf &&
    node->parent->SEChild->isLeaf;
}

bool Neighbor::levelPassed(Node * node){
    vector<Node*> vec;
    vector<vector<Node*> > neighbors (4, vec);
    getNeighbors(node,neighbors);
    int newLevel = node->currentLevel+1;
    for(vector<vector<Node*> >::iterator IT = neighbors.begin();
        IT!=neighbors.end();IT++){
        for(vector<Node*>::iterator direction = (*IT).begin();
            direction!=(*IT).end();direction++){
            if(abs(newLevel - (*direction)->currentLevel)>1)
                return false;
        }
    }
    return true;
    
}