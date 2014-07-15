//
//  See treeRenderer.h for description
//
//
//  Created by Alexandra Gendreau on 1/21/14.
//
//


#include "treeRenderer.h"

using namespace std;

/*
 * Render text to screen, method found via Google search
 */

 void treeRenderer::drawString (void * font, std::string str,
                               double x, double y, double z){
    unsigned int i;
    glRasterPos3f(x, y, z);
    
    for (string::iterator it=str.begin(); it!=str.end(); it++)
        glutBitmapCharacter (font, *it);
}


/* 
 * Draw filled rectangles that represent the neighbors of the given node
 * currently colors preset, this could change
 * red = north, cyan = south, blue = east, magenta = west
 * We could use the rectangle class mentioned below here as well
 */

void treeRenderer::displayHelperNeighbors(std::vector<std::vector<
                                          Rectangle> > neighborInfo){
    double x,y,w,h;
    int i   = 0; //for coloring
    
    for(vector<vector<Rectangle> >::iterator directionIterator = neighborInfo.begin();
        directionIterator!=neighborInfo.end();++directionIterator){
        glColor3f (1.0, 0.0, 0.0);
        switch(i) {//set colors neighbors are stored north, south, east, west
            case 0:
                glColor3f (1.0, 0.0, 0.0);
                break;
            case 1:
                glColor3f (0.0, 1.0, 1.0);
                break;
            case 2:
                glColor3f (0.0, 0.0, 1.0);
                break;
            case 3:
                glColor3f (1.0, 0.0, 1.0);
                break;
        }
        for(vector<Rectangle>::iterator neighborIterator = (*directionIterator).begin();
            neighborIterator!= (*directionIterator).end();++neighborIterator){
            x   = neighborIterator->leftX; //lower left x
            y   = neighborIterator->leftY; //lower left y
            w   = neighborIterator->width; //width
            h   = neighborIterator->height; //height
            glRectf(x,y,x+w,y+h);
        }
        i++;
    }
}


/* 
 * Draw the spatial division represented by the quad tree.
 * a vector of vectors is passed
 *
 * Wach vector holds the lower left coordinates of a rectangle as well as
 * its height and width
 * Using a vector of vectors might not be the best option.  Maybe I should create
 * a little rectangle class or something
 */


void treeRenderer::displayHelperTree(std::vector<Rectangle> gridInfo){
    
    double x,y,w,h;
    glColor3f (0.0, 1.0, 0.0);
    
    for(vector<Rectangle>::iterator directionIterator = gridInfo.begin();
        directionIterator!=gridInfo.end();directionIterator++){
        x   = directionIterator->leftX; //lower left x
        y   = directionIterator->leftY; //lower left y
        w   = directionIterator->width; //width
        h   = directionIterator->height; //height
        
        glBegin(GL_LINE_LOOP);
        glVertex2f (x, y);
        glVertex2f (x+w, y);
        glVertex2f (x+w, y+h);
        glVertex2f (x, y+h);
        glEnd();
    }
    
}


/* 
 * draws a straight line from (x0,y0) to (x1,y1), currently colore black
 */

 void treeRenderer::drawLine(double x0,double y0,double x1,double y1){
    
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINES);
    glVertex2f(x0,y0);
    glVertex2f(x1,y1);
    glEnd();
    
}
