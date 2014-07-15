//
//  OneLevel.cpp
//  
//
//  Created by Alexandra Gendreau on 2/13/14.
//
//

#include "OneLevel.h"

OneLevel::OneLevel(double x,double y, double width, double height,
                   int numCells,int max,Application * app) : QuadTree(x,y,width,height, numCells, max, app){
}

OneLevel::~OneLevel(){}

bool OneLevel::refine(Node * node){
    return QuadTree::refine(node);
}

bool OneLevel::coarsen(Node * node){
    
    return QuadTree::coarsen(node) &&
    node->parent->NEChild->isLeaf &&
    node->parent->NWChild->isLeaf &&
    node->parent->SWChild->isLeaf &&
    node->parent->SEChild->isLeaf;
}