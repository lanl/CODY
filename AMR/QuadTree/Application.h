//
//  Application.h
//  The application determines the local refinement and coarsening criteria
//
//  Abstract base the class which allows users to implemented their own subclass
//  with their given refinement/coarsening criteria
//
//  Created by Alexandra Gendreau on 2/11/14.
//
//
//

#ifndef _Application_h
#define _Application_h

#include <iostream>
#include <memory>

#endif

class Application {
public:
 
    virtual ~Application(){}
    
    virtual bool refine(double x, double y, double w, double h) = 0;
    virtual bool coarsen(double x, double y, double w, double h) = 0;
    

};

/*
 * This class is used for interacting with the quadtree visualization
 * where clicking on a node always coarsens or refines
 */

class Interaction : public Application {
public:
    
    Interaction(){
    }
    ~Interaction(){}
    bool refine(double x, double y, double w, double h){
        return true;
    }
    
    bool coarsen(double x, double y, double w, double h){
        return true;
    }
    
};


/*
 * Line class
 * Intersection with the line is the refinement criteria for an application where
 * a line is "pushed" through the grid
 * We want to line to extend across the entire
 * space, xS = the left x boundary and xE = the right x boundary
 */

class Segment { 
    
public:
    
    Segment(double xS,double xE,double m,double b){
        //start coordinates
        x0          = xS;  
        y0          = m*x0+b;
        
        //finish coordinates
        x1          = xE;
        y1          = m*x1+b;
        
        slope       = m;
        intercept   = b;
    }
    
    /*
     * Translates the line x in the horizontal direction and 
     * y in the vertical direction
     */
    void translate(double x,double y){
        intercept   = intercept - (slope*x)+y;
        y0          = slope*x0+intercept;
        y1          = slope*x1+intercept;
    }
    
    /*
     * Resets the line back to the original
     */
    void reset(double originalB){
        intercept   = originalB;
        y0          = slope*x0+intercept;
        y1          = slope*x1+intercept;
    }
    
    void setx0(double x){
        x0  = x;
        y0  = slope*x0+intercept;
    }
    
    
    /*
     * Getters and Settings
     */
    void setx1(double x){
        x1  = x;
        y1  = slope*x1+intercept;
    }
    
    double getx0(){
        return x0;
    }
    
    double getx1(){
        return x1;
    }
    
    double gety0(){
        return y0;
    }
    
    double gety1(){
        return y1;
    }
    
private:
    /*
     * Instance variables 
     */
    double x0,y0,x1,y1,slope,intercept;
    
};

/*
 * This class is used for pushing a line through a grid and refining when the 
 * line and the node intersect 
 */

class Line : public Application {
public:
    Line(Segment * segment){
        seg     = segment;
    }
    
    ~Line(){
    }
    
    /*
     * Intersection test
     */
    bool refine(double x, double y, double w, double h){
        double minX     = seg->getx0();
        double maxX     = seg->getx1();
        double dx       = seg->getx1() - seg->getx0();;
        
        double minY     = seg->gety0();
        double maxY     = seg->gety1();
        
        //x projection
        if(maxX > (x+w))
            maxX = x+w;
        if(minX < x)
            minX = x;
        if(minX >= maxX)
            return false;
        
        //y projection
        if(dx > 0){
            double a        = (seg->gety1()-seg->gety0())/dx;
            double b        = seg->gety0() - a*seg->getx0();
            minY            = a*minX +b;
            maxY            = a*maxX +b;
            if(minY > maxY){
                double temp = maxY;
                maxY        = minY;
                minY        = temp;
            }
        }
        if(maxY > (y+h))
            maxY = y+h;
        if(minY < y)
            minY = y;
        if(minY >= maxY)
            return false;
        
        return true;//the segment intersects the node, so we refine if possible
    }
    
    bool coarsen(double x, double y, double w, double h){
        return !(refine(x,y,w,h)); //the segment does not intersect the node
    }
    
    Segment * getSegment(){
        return seg;
    }
    
private:
    Segment * seg;
};
