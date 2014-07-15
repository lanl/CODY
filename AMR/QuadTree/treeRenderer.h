//
//  This class contains the methods that draw the quad tree and text
//  I made these static methods because I could not think of any necessary
//  instance variables, so I thought it made more sense to just make them class
//  methods.  If at same point I give treeRenderer objects instance variables,
//  I will probably change these methods to instance methods.
//
//  Created by Alexandra Gendreau on 1/21/14.
//
//

#ifndef ____treeRenderer__
#define ____treeRenderer__

#include <cstdio>
#include "GLUT/glut.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <iostream>

/* 
 * Simple struct to hold coordinates of quad tree leaves
 */
struct Rectangle{
    double leftX,leftY, width, height;
    Rectangle(double x, double y, double w, double h){
        leftX   = x;
        leftY   = y;
        width   = w;
        height  = h;
    }
};

class treeRenderer{
public:
    
    /*
     * Render text to screen, method found via Google search
     */
     void static drawString (void * font, std::string str,
                            double x, double y, double z);
    
    /* 
     * Draw filled rectangles that represent the neighbors of the given node
     * currently colors preset, this could change
     * red = north, cyan = south, blue = east, magenta = west
     * We could use the rectangle class mentioned below here as well
     */
    void static displayHelperNeighbors(std::vector<std::vector
                                       <Rectangle> > neighborInfo);
    
    /* 
     * Draw the spatial division represented by the quad tree.
     * a vector of vectors is passed
     *
     * Each vector holds the lower left coordinates of a rectangle as well as
     * its height and width
     * Using a vector of vectors might not be the best option.  Maybe I should create
     * a little rectangle class or something
     */
    void static displayHelperTree(std::vector<Rectangle> gridInfo);
    
    /*
     * Draws a straight line from (x0,y0) to (x1,y1)
     */
    void static drawLine(double x0, double y0, double x1, double y1);
    
};

#endif /* defined(____treeRenderer__) */