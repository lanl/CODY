import std.stdio;
import std.parallelism;
import core.thread;
import std.string;
import std.datetime;
import std.getopt;
import std.conv;
import std.math;

//Global Variables

TaskPool	workers;   //global variable for workers, allow overloading
shared Segment 	seg; 	   //shared segment, should each node hold it's own?


/*
 * This is the basic class used to build the quad tree
 * It contains pointers to each of its children (null if leaf)
 * A pointer to its parent
 * Which child of its parent (root defaults to -1)
 * A Geom to specify its location in space
 * Its level in the tree
 * A dynamic array for its north, south, east, and west neighbors
 * 
 * A bool to determine if we should coarsen the node
 * used during concurrent update
 */


class Node { 
	Node 		NEChild;
	Node		NWChild;
	Node		SWChild;
        Node		SEChild;
	Node 		parent;
	Geom 		geom;
	int 		currentLevel;
	int 		childType;
	bool 		isLeaf;
	bool 		toCoarsen;
	Node[][] 	neighbors; 	//holds neighbors
	
	//constructor
	this(Geom geo,Node par,int level,int cType) {
		parent		= par;
		geom 		= geo;
		currentLevel 	= level;
		childType 	= cType;
		isLeaf 		= true;
		toCoarsen 	= false;

		Node[][] n;
		for(int i=0;i<4;i++){
			Node[] neigh;
			n=n~neigh;
		}

		neighbors 	= n;
	}
}


/*
 * This struct holds the geometric information for each node of the tree
 * Each node is a rectangle where
 * x,y is the lower left corner
 * w,h are the width and height respectively
 */

struct Geom {

	double	x,y;
	double	w,h;

	//Constructor
	this(double xStart, double yStart, double width, double height) {
		x 	= xStart;
		y	= yStart;
		w 	= width;
		h 	= height;
	}

	//method for printing the geometry
	void printGeom() {
		writeln("(",x,",",y,") ");
		writeln("height: ",h," width: ",w);
	}
}


/*
 * This struct holds the informatiob for the line used in the refinement and
 * coarsening process.  The line consists of 
 * x0,y0 is the initial point
 * x1,y1 is the final point
 * slope is the slope of the line
 * intercept is the y-intercept of the line
 */

struct Segment {
	double 	x0,y0;
	double 	x1,y1;
	double 	slope, intercept;

	//Constructor
	this(double xStart,double xEnd, double m, double b){
		x0 		= xStart;
		x1 		= xEnd;
		slope 		= m;
		intercept 	= b;

		y0 		= x0*slope+intercept;
		y1 		= x1*slope+intercept;
	}
	
	/* 
	 * Method that translates a line segment by
	 * d1 in the x-direction
	 * d2 in the y-direction
	 * shared in order to be used concurrently
	 */
	shared void translate(double d1,double d2) {
		intercept 	= intercept - (slope * d1) + d2;
		y0 		= slope*x0 + intercept;
		y1 		= slope*x1 + intercept;
	}
}

class QuadTree {//Quad tree
	Node 	root;		//reference to root
				//reference since Node is a class
	
	int 	maxLevel;	//maximum number of levels in tree
				//root is level 1

	//Constructor
	this(double x, double y, double width, double height,int levels,
			int max) {//,Segment seg) {
		maxLevel	= max;
		Geom geom 	= Geom(x,y,width,height);
		root 		= insert(null,geom,1,-1,levels);//,seg);
	}

	/* 
	 * This method constructrs the complete quad tree to the deseired level
	 * The pointer to the root is returned
	 */
	Node insert(Node par,Geom geom,int level, int CT, 
			int levelsRemaining) {//,Segment seg) {
		if(levelsRemaining == 1){ //reached the leaf level
			Node node 	= new Node(geom,par,level,CT);//,seg);
		       	return node;

		} else { //internal node
			Node node 	= new Node(geom,par,level,CT);//,seg);
			node.isLeaf 	= false;

			//recursively construct the NE child
			Geom newGeom 	= Geom(geom.x + (geom.w/2.0),
						geom.y + (geom.h/2.0),geom.w/2.0,
						geom.h/2.0);
			node.NEChild 	= insert(node,newGeom,level + 1,0,
						levelsRemaining - 1);//,seg);

			//recursively construct the NW child
			newGeom.x 	= geom.x;
			newGeom.y 	= geom.y + newGeom.h;
			node.NWChild 	= insert(node,newGeom,level + 1,1,
						levelsRemaining - 1);//,seg);

			//recursively construct the SW child
			newGeom.x 	= geom.x;
			newGeom.y 	= geom.y;
			node.SWChild 	= insert(node,newGeom,level + 1,2,
						levelsRemaining - 1);//,seg);
			
			//recursively construct the SE child
			newGeom.x 	= geom.x + newGeom.w;
			newGeom.y 	= geom.y;
			node.SEChild 	= insert(node,newGeom,level + 1,3,
						levelsRemaining - 1);//,seg);
			
			return node;
		}
	}

	
	/*
	 * This method determines if the segment and node intersect
	 * Note it does not return where they intersect just if they intersect
	 */
	
	bool intersect(Node node,Segment seg) {
		double minX 	= seg.x0;
		double maxX 	= seg.x1;
		double dx 	= maxX - minX;

		double minY 	= seg.y0;
		double maxY 	= seg.y1;

		if(maxX	> (node.geom.x + node.geom.w)) {
			maxX 	= node.geom.x + node.geom.w;
		}

		if(minX < node.geom.x) {
			minX 	= node.geom.x;
		}

		if(minX >= maxX) {
			return false;
		}

		if(dx > 0){
			minY 	= seg.slope * minX + seg.intercept;
			maxY 	= seg.slope * maxX + seg.intercept;

			if(minY > maxY) {
				auto temp 	= maxY;
				maxY 		= minY;
				minY 		= temp;
			}
		}

		if(maxY > (node.geom.y + node.geom.h)) {
			maxY = node.geom.y + node.geom.h;
		}
		if(minY < node.geom.y){
			minY = node.geom.y;
		}
		if(minY >= maxY) {
			return false;
		}

		return true;
	}

	/*
	 * This method refines the node using the line segment criteria
	 * Recursively refines as deep as possible
	 */

	void refineNode(Node node) {
		Geom geom 	= node.geom;
		int level 	= node.currentLevel + 1;
		node.isLeaf 	= false;
		//Segment seg 	= node.segment;


		//construct NE child
		Geom newGeom 	= Geom(geom.x + (geom.w/2.0),
					geom.y + (geom.h/2.0),
					geom.w/2.0,geom.h/2.0);
		node.NEChild 	= new Node(newGeom,node,level,0);//,seg);

		//construct NW child
		newGeom.x 	= geom.x;
		newGeom.y 	= geom.y + newGeom.h;
		node.NWChild 	= new Node(newGeom,node,level,1);//,seg);
		

		//construct SW child
		newGeom.x 	= geom.x;
		newGeom.y 	= geom.y;
		node.SWChild 	= new Node(newGeom,node,level,2);//,seg);

		//construct SE child
		newGeom.x 	= geom.x + newGeom.w;
		newGeom.y 	= geom.y;
		node.SEChild 	= new Node(newGeom,node,level,3);//,seg);

		//Recursively keep refining
		if(refine(node.NEChild)) {
			refineNode(node.NEChild);
		}

		if(refine(node.NWChild)) {
			refineNode(node.NWChild);
		}

		if(refine(node.SWChild)) {
			refineNode(node.SWChild);
		}

		if(refine(node.SEChild)) {
			refineNode(node.SEChild);
		}
	}

	
	/*
	 * This method determines whether or not to refine a node
	 * Refines if the level is less than the maximum and
	 *	the node and segment intersect
	 */

	bool refine(Node node) {
		bool value 	= node.currentLevel < maxLevel &&
					intersect(node,seg);
		return value;
	}

	
	/*
	 * This method coarsens the node by setting its children to null
	 *	The children are garbage collected
	 *	The node is now tagged as a leaf node
	 */

	void coarsenNode(Node node) {
		node.NEChild 	= null;
		node.NWChild 	= null;
		node.SWChild 	= null;
		node.SEChild 	= null;
		node.isLeaf 	= true;
		node.toCoarsen 	= false;
	}

	
	/*
	 * This method determines whether or not to coarsen a node
	 * Coarsens if the line segment does not intersect the node
	 * Only called on nodes that are not leaves
	 */

	bool coarsen(Node node) {
		return !(intersect(node,seg));
	}

		
	/*
	 * This method calls the other methods to the fill the neighbor 
	 * 	dynamic arrays for each node
	 */

	void getNeighbors(Node node, ref Node[][] neighbors){
		
		if(node.parent is null) { //parents have no neighbors
			return;
		}

		//Get the neighbor array of the parent node
		Node[][] 	parentsNeighbors;
	        for(int i=0;i<4;i++){
			Node[] 			neigh;
	       		parentsNeighbors 	= parentsNeighbors~neigh;
		}

		getParentsNeighbors(parentsNeighbors,node.parent);

		//How the neighbor array is filled depends on its child type
		switch(node.childType) {
			case 0: //NE Child
				getNeighborsSibs(node.parent.NWChild,
						neighbors[3],0,3);
				getNeighborsSibs(node.parent.SEChild,
						neighbors[1],0,1);
				processNeighbors(parentsNeighbors[2],
						neighbors[2],
						node.currentLevel-1,1);
				processNeighbors(parentsNeighbors[0],
						neighbors[0],
						node.currentLevel-1,3);
				break;

			case 1: //NW Child
				getNeighborsSibs(node.parent.NEChild,
						neighbors[2],1,2);
				getNeighborsSibs(node.parent.SWChild,
						neighbors[1],0,1);
				processNeighbors(parentsNeighbors[3],
						neighbors[3],
						node.currentLevel-1,0);
				processNeighbors(parentsNeighbors[0],
						neighbors[0],
						node.currentLevel-1,2);
				break;

			case 2: //SW Child
				getNeighborsSibs(node.parent.SEChild,
						neighbors[2],1,2);
				getNeighborsSibs(node.parent.NWChild,
						neighbors[0],2,3);
				processNeighbors(parentsNeighbors[3],
						neighbors[3],
						node.currentLevel-1,3);
				processNeighbors(parentsNeighbors[1],
						neighbors[1],
						node.currentLevel-1,1);
				break;

			case 3: //SE Child
				getNeighborsSibs(node.parent.SWChild,
						neighbors[3],0,3);
				getNeighborsSibs(node.parent.NEChild,
						neighbors[0],2,3);
				processNeighbors(parentsNeighbors[2],
						neighbors[2],
						node.currentLevel-1,2);
				processNeighbors(parentsNeighbors[1],
						neighbors[1],
						node.currentLevel-1,0);
				break;

			default://d forces a default case
				//this should never match
				writeln("SOMETHING REALLY BAD");
				break;
		}
	}
	

	/*
	 * This function fills the dynamic array of the parent's neighbors
	 */
	
	void getParentsNeighbors(Node[][] parentsNeighbors,Node node) {
		getNeighbors(node, parentsNeighbors);
	}

	
	/*
	 * This function returns the neighbors that are either the sibling
	 * of the node or the children of a sibling since all neighbors
	 * are leaf nodes
	 * 
	 * Based on the child type of the siblings, the function recurses
	 * down the tree until it reaches the correct neighbor leaves 
	 * and adds them to the dynamic array
	 */

	void getNeighborsSibs(Node node,ref Node[] list,
			int sibOneType,int sibTwoType) {
		if(node.isLeaf) {
			list 	= list~node;
		} else {
			if(sibOneType == 0 || sibTwoType == 0) {
				getNeighborsSibs(node.NEChild,list,
						sibOneType,sibTwoType);
			}

			if(sibOneType == 1 || sibTwoType == 1) {
				getNeighborsSibs(node.NWChild,list,
						sibOneType,sibTwoType);
			}

			if(sibOneType == 2 || sibTwoType == 2) {
				getNeighborsSibs(node.SWChild,list,
						sibOneType,sibTwoType);
			}
			if(sibOneType == 3 || sibTwoType == 3) {
				getNeighborsSibs(node.SEChild,list,
						sibOneType,sibTwoType);
			}
		}
	}


	/*
	 * This function determines whether or not the node is actually a 
	 * neighbor
	 * 
	 * Based on the processing of the neighbors of the parent
	 */

	bool actualNeighbors(Node node,int level, int nodeType) {
		if(node.currentLevel == level) {
			return node.childType == nodeType;
		} else {
			return actualNeighbors(node.parent,level,nodeType);
		}
	}


	/*
	 * This function processes the neighbor array of the parent node 
	 * to determines which of thos nodes are also neighbors
	 */

	void processNeighbors(Node[] pNeighbors,ref Node[] list,int pLevel,
			int nodeType) {
		if(pNeighbors.length == 0) {
			return;
		} else if(pNeighbors[0].currentLevel<=pLevel) {
			list 	= list~pNeighbors[0];
		} else {
			for(int i = 0;i<pNeighbors.length;i++) {
				if((pNeighbors[i].currentLevel == (pLevel+1)) &&
						(pNeighbors[i].childType ==
						 nodeType)) {
					list 	= list~pNeighbors[i];
				} else
					if(actualNeighbors(pNeighbors[i],pLevel+1,nodeType)){
						list 	= list~pNeighbors[i];
					}
			}
		}
	}


/************************** TESTING AND UPDATING ***********************************/

	/*
	 * This method updates the quad tree serially
	 * Used in timing tests 
	 */
	
	void updateSerial(Node node) {
		if(node.isLeaf) {
			if(refine(node)) {
				refineNode(node);
			}

		} else if(coarsen(node)){
			coarsenNode(node);

		} else {
			updateSerial(node.NEChild);
			updateSerial(node.NWChild);
			updateSerial(node.SWChild);
			updateSerial(node.SEChild);
		}
	}

	
	/*
	 * This method is called by each task when updating the quad tree
	 * concurrently
	 */

	void update(Node node) {
		if(refine(node)) {
			refineNode(node);
		} else if(node.toCoarsen) {
			coarsenNode(node);
		}
	}

	
	/*
	 * This method adds work to the task pool for updating
	 */

	void work(Node node, void delegate(Node) func) {
		if(node.isLeaf) {	
			auto t 	= task(func,node);
			workers.put(t);
		} else if(coarsen(node)){//we should coarsen
			//prevents continued traversal
			node.toCoarsen 	= true;
			auto t 		= task(func,node);
			workers.put(t);
		} else {
			work(node.NEChild,func);
			work(node.NWChild,func);
			work(node.SWChild,func);
			work(node.SEChild,func);
		}
	}

	/*
	 * Adds work to the taskpool to concurrently find neighbors of each node
	 */
	
	void neighborWork(Node node, void delegate(Node, ref Node[][]) func) {
		if(node.isLeaf) {
			auto t 	= task(func,node,node.neighbors);
			workers.put(t);
		} else {
			auto t 	= task(func,node,node.neighbors);
			workers.put(t);
			neighborWork(node.NEChild,func);
			neighborWork(node.NWChild,func);
			neighborWork(node.SWChild,func);
			neighborWork(node.SEChild,func);
		}
	}

	
	/*
	 * This method executes dummy work on each leaf node
	 * Helps to determine how much work is need to justify concurrency
	 */

	void dumby(Node node, int iter) {
		for(int i = 0; i < iter; i++) {
			//dummy work
		}
	}


	/*
	 * This executes the dummy work serially
	 * Used for timing test
	 */

	void dumbyWorkSerial(Node node,int iter,void delegate(Node,int) func){
		if(node.isLeaf) {
			func(node,iter);
		} else {
			dumbyWorkSerial(node.NEChild,iter,func);
			dumbyWorkSerial(node.NWChild,iter,func);
			dumbyWorkSerial(node.SWChild,iter,func);
			dumbyWorkSerial(node.SEChild,iter,func);
		}
	}


	/*
	 * This method adds dummy work for each leaf node to the task pool
	 * Used in timing tests
	 */
	
	void dumbyWork(Node node, int iter,void delegate(Node,int) func) {
		if(node.isLeaf) {	
			auto t 	= task(func,node,iter);
			workers.put(t);
		} else {
			dumbyWork(node.NEChild,iter,func);
			dumbyWork(node.NWChild,iter,func);
			dumbyWork(node.SWChild,iter,func);
			dumbyWork(node.SEChild,iter,func);
		}
	}


/************************ DEBUGGING STUFF ****************************/

	/*
	 * This method returns how many nodes are in the tree
	 * Invokes the helper method countNodesHelper
	 */

	int countNodes() {
		return countNodesHelper(root);
	}

	
	/*
	 * This method is a helper method for counting the nodes in the tree
	 */

	int countNodesHelper(Node node) {
		if(node.isLeaf) {
			return 1;
		} else {
			return 1 +
				countNodesHelper(node.NEChild) +
				countNodesHelper(node.NWChild) +
				countNodesHelper(node.SWChild) +
				countNodesHelper(node.SEChild);
		}
	}

	
	/*
	 * This method prints the level of each 
	 * Invokes the helper method printLevelsHelper
	 */

	void printLevels() {
		return printLevelsHelper(root);
	}


	/*
	 * This method prints the levels of each node
	 * Pre-Order 
	 */
	void printLevelsHelper(Node node) {
		if(node.isLeaf) {
			writeln(node.currentLevel);
		} else {
			writeln(node.currentLevel);
			printLevelsHelper(node.NEChild);
			printLevelsHelper(node.NWChild);
			printLevelsHelper(node.SWChild);
			printLevelsHelper(node.SEChild);
		}
	}

	
	/*
	 * This method prints the geometry of each node
	 * Invokes the helper method printGeomHelper
	 */
	void printGeom() {
		return printGeomHelper(root);
	}


	/*
	 * This method prints the geomtry of each node
	 * Uses the printGeom method of the Geom struct
	 */

	void printGeomHelper(Node node) {
		if(node.isLeaf) {
			node.geom.printGeom();
		} else {
			node.geom.printGeom();
			printGeomHelper(node.NEChild);
			printGeomHelper(node.NWChild);
			printGeomHelper(node.SWChild);
			printGeomHelper(node.SEChild);
		}
	}

	
	/*
	 * This method counts the number of leaves in the tree
	 * Invokes the helper method countLeavesHelper
	 */
	
	int countLeaves() {
		return countLeavesHelper(root);
	}


	/*
	 * This method counts the number of leaves in the tree
	 */

	int countLeavesHelper(Node node) {
		if(node.isLeaf) {
			return 1;
		} else {
			return countLeavesHelper(node.NEChild)+
				countLeavesHelper(node.NWChild)+
				countLeavesHelper(node.SWChild)+
				countLeavesHelper(node.SEChild);
		}
	}

	
	/*
	 * This method outputs the area of the leaves
	 * Used for checking that the tree was updated correctly
	 * Invokes the helper method areaTestHelper
	 */

	float areaTest() {
		return areaTestHelper(root);
	}


	/*
	 * This method returns the area of the leaves
	 */

	float areaTestHelper(Node node){
		if(node.isLeaf) {
			return node.geom.w * node.geom.h;
		} else {
			return areaTestHelper(node.NEChild)+
				areaTestHelper(node.NWChild)+
				areaTestHelper(node.SWChild)+
				areaTestHelper(node.SEChild);
		}
	}

	
	
	}

void main(string[] args) {
	
	QuadTree test = new QuadTree(-4.0,-4.0,4.0,4.0,3,3);//,seg);
	
	//initialize the segment and stopwatch for timing
	double originalB = -7.0;
	seg = Segment(-4.0,0.0,-1.0,originalB);
	seg.translate(1.5,1.5);
	StopWatch sw;


	//set up command line arguments
	int funcType;
	int depth; 
	int maxIter;
	int dumbyIter; 
	string filename;
	int numCores;

	//initialize command line arguments
	getopt(args,
			"case", &funcType,
			"filename", &filename,
			"depth", &depth,
			"maxIter", &maxIter,
			"dumbyIter", &dumbyIter,
			"numCores", &numCores);


	//run tests
	switch(funcType){
		case 0: //concurrent work test
			writeln("Concurrent Test\n");
			auto f 	= File(filename,"w");
			f.write("leaves,max depth,time,traversal,nodes\n");
			for(int i = 8; i < depth; i=i+2) {
				//build and update tree, check correctness
				QuadTree qTree= new 
					QuadTree(-4.0,-4.0,4.0,4.0,i,i);
				/*qTree.updateSerial(qTree.root);
				if(qTree.areaTest() != (qTree.root.geom.w *
						 qTree.root.geom.h)){
					writeln("Actual Area: ",qTree.areaTest());
					break;
				}*/
				StopWatch w;	
				w.start();
				for(int j = 0; j < maxIter; j++) {
					workers 	= new TaskPool(numCores-1);
					qTree.dumbyWork(qTree.root,
							dumbyIter,&qTree.dumby);
					
					workers.finish(true);
				}
				TickDuration time 	= w.peek();
				//sw.reset();
				w.stop();
				//time traversal with no work

				StopWatch trav;
				trav.start();
				for(int k = 0; k<maxIter; k++){
					workers 	= new TaskPool(numCores-1);
					qTree.dumbyWork(qTree.root,1,&qTree.dumby);
					workers.finish(true);
				}
				TickDuration tTime 	= trav.peek();
				//sw.reset();
				trav.stop();
				
				
				//writeln("traversal time: ", tTime.usecs/to!float(maxIter));
				//time the dummy work on leaves
				//ouput results to file
				f.write(pow(4,i-1),",",i,",",
						(time.usecs)/to!float(maxIter),",",
						(tTime.usecs)/to!float(maxIter),",",
						(pow(4,i)-1)/3,"\n");
			}
			break;

		case 1: //sequential depth test
			writeln("Serial Test\n");
			auto f 	= File(filename,"w");
			f.write("leaves,max depth,time,traversal,nodes\n");
			for(int i = 8; i < depth; i=i+2){
				//build and update tree, check correctness
				QuadTree qTree = new
					QuadTree(-4.0,-4.0,4.0,4.0,i,i);
				/*qTree.updateSerial(qTree.root);
				if(qTree.areaTest() != (qTree.root.geom.w *
							qTree.root.geom.h)){
					writeln("Actual Area: ",qTree.areaTest());
					break;
				}*/
				StopWatch w;
				w.start();
				for(int j = 0 ;j < maxIter; j++){
					qTree.dumbyWorkSerial(qTree.root,
							dumbyIter,&qTree.dumby);
				}
				TickDuration time 	= w.peek();
				w.stop();
				//time traversal with no work
				//sw.reset();
				
				StopWatch trav;
				trav.start();
				for(int k = 0; k<maxIter; k++){
					qTree.dumbyWorkSerial(qTree.root,1,&qTree.dumby);
				}

				TickDuration tTime 	= trav.peek();
				//sw.reset();
				trav.stop();

				
				//writeln("traversal time: ", tTime.usecs/to!float(maxIter));
				//time the sequential dummy work on leaves
				//output results to file
				f.write(pow(4,i-1),",",i,",",
						(time.usecs)/to!float(maxIter),",",
						(tTime.usecs)/to!float(maxIter),",",
						(pow(4,i)-1)/3,"\n");
			}
			break;
		case 2:
			int work = dumbyIter;
			//time for dummy work
			sw.start();
			for(int t=0;t<1000;t++){
			test.dumby(test.root,work);
			}
			sw.stop();
			TickDuration timeDumby = sw.peek();
			writeln("Iterations: ",work);
			writeln("Time: ",timeDumby.usecs/to!float(1000));
			sw.reset();
			break;
		default:
			writeln("OOPS, bad function case");
			break;		
	}
}
