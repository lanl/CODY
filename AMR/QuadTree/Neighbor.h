//
//  Neighbor.h
//  
//
//  Subclass of QuadTree which will not allow a nodes to have more than one
//  level different of refinement.
//
//  Achieved by overriding the global refine/coarsen methods in the QuadTree
//  base class
//
//

#ifndef ____Neighbor__
#define ____Neighbor__

#include <iostream>
#include "QuadTree.h"

class Neighbor : public QuadTree {
    
public:
    Neighbor(double x,double y, double width, double height,
             int numCells,int max,Application * app);
    ~Neighbor();
    
    bool refine(Node * node);
    bool coarsen(Node * node);
    
private:
    bool levelPassed(Node * node);
};

#endif /* defined(____Neighbor__) */