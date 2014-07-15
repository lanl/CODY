//
//  OneLevel.h
//  
//
//  Subclass of QuadTree which will only allow nodes to be coarsened if
//  all of their children are leaves.
//
//  Used for visualizations
//
//  Achieved by overriding the global refine/coarsen methods in the QuadTree
//  base class
//
//

#ifndef ____OneLevel__
#define ____OneLevel__

#include <iostream>
#include "QuadTree.h"

class OneLevel : public QuadTree {
    
public:
    OneLevel(double x,double y, double width, double height,
             int numCells,int max,Application * app);
    ~OneLevel();
    
    virtual bool refine(Node * node);
    virtual bool coarsen(Node * node);
};

#endif /* defined(____OneLevel__) */