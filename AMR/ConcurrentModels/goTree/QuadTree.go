package main

import (
	"fmt"
	"runtime"
	"time"
	"os"
	"flag"
)


//Global Variables

var maxLevel 	int		//maximum number of levels in tree
var root	*Node		//pointer to root of tree
var segment 	*Segment	//segment used for refinement/coarsening


/*
 * This struct hold the geometric information for each node of the tree
 * Each node is a rectangle where
 * x,y is the lower left corner
 * w,h are the width and height respectively 
 */

type Geom struct {
	x,y	float32
	w,h	float32
}


/*
 * This struct hold the information for the line used in the refinement
 * and coarsening process.  The line consists of
 * x0,y0 is the initial point 
 * x1,y1 is the final point
 * slope is the slope of the line
 * intercept is the y-intercept of the line
 */

type Segment struct {
	x0,y0			float32
	x1,y1	  		float32
	slope,intercept 	float32
}


/*
 * Method that output the segment's initial and final points to a file
 */

func (seg *Segment) outputSegment(file *os.File) {
	line := fmt.Sprintf("%f %f \n %f %f \n",seg.x0,seg.y0,seg.x1,seg.y1)
	l    := []byte(line);
	file.Write(l);
}


/*
 * Method that translate's a line segment by 
 * f1 in the x-direction
 * f2 in the y-direction
 */

func (seg *Segment) Translate(f1,f2 float32){
	seg.intercept = seg.intercept - (seg.slope * f1) + f2
	seg.y0	      = seg.slope*seg.x0 + seg.intercept
	seg.y1 	      = seg.slope*seg.x1 + seg.intercept
}


/*
 * This is the basic struct used to build the quad tree
 * It contains pointers to each of its children (nil if leaf)
 * A pointer to its parent
 * Which child of its parent (root defaults to -1)
 * A Geom to specify its location in space
 * Its level in the tree
 * A slice of its north, south, east, and west neighbors
 * 
 * A bool to determine if we should coarsen the node
 * used during concurrent update
 */

type Node struct {
	NE 		*Node
	NW 		*Node
	SW 		*Node
	SE 		*Node
	parent 		*Node
	isLeaf 		bool
	geo 		Geom
	currentLevel 	int
	childType 	int
	neighbors 	[][]*Node
	toCoarsen 	bool
}


/*
 * This function constructs the complete  quad tree to the desired level
 * The pointer to the root is returned
 */

func Construct(par *Node, geom Geom,level,CT,levelRemaining int) *Node {
	if levelRemaining == 1 { //reached the leaf level
		node			:= new (Node)
		node.geo 		= geom
		node.parent 		= par
		node.isLeaf 		= true
		node.currentLevel 	= level
		node.childType 		= CT
		return node
	} else { //internal node
		node			:=new (Node)	
		node.parent 		= par
		node.currentLevel 	= level
		node.childType 		= CT
		node.isLeaf 		= false
		node.geo 		= geom

		//recursively construct NE child
		newGeom			:= Geom{}
		newGeom.x 		= geom.x+geom.w/2.0
		newGeom.y 		= geom.y+geom.h/2.0
		newGeom.w 		= geom.w/2.0
		newGeom.h 		= geom.h/2.0
		node.NE 		= Construct(node,newGeom,level+1,0,levelRemaining-1)

		//recursively construct NW child
		newGeom.x 		= geom.x
		newGeom.y 		= geom.y + newGeom.h
		node.NW 		= Construct(node,newGeom,level+1,1,levelRemaining-1)
		
		//recursively construct SW child
		newGeom.x 		= geom.x
		newGeom.y 		= geom.y
		node.SW 		= Construct(node,newGeom,level+1,2,levelRemaining-1)
		
		//recursively construct SE child
		newGeom.x 		= geom.x+newGeom.w
		newGeom.y 		= geom.y
		node.SE 		= Construct(node,newGeom,level+1,3,levelRemaining-1)

		return node
	}
}


/*
 * This function refines the node using the line segment criteria
 * Recursively refines as deep as possible
 */

func refineNode(node *Node) {
	//update internal node	
	geom 		:= node.geo
	node.isLeaf 	= false	
		
	//construct NE child
	newGeom		:= Geom{}
	newGeom.x 	= geom.x+geom.w/2.0
	newGeom.y 	= geom.y+geom.h/2.0
	newGeom.w 	= geom.w/2.0
	newGeom.h 	= geom.h/2.0
	node.NE 	= new(Node)
	newNode(node.NE,node,newGeom,0)

	//construct NW child
	newGeom.x 	= geom.x
	newGeom.y 	= geom.y+newGeom.h
	node.NW 	= new(Node)
	newNode(node.NW,node,newGeom,1)

	//construct SW child
	newGeom.x 	= geom.x
	newGeom.y 	= geom.y
	node.SW 	= new(Node)
	newNode(node.SW,node,newGeom,2)

	//construct SE child
	newGeom.x 	= geom.x+newGeom.w
	newGeom.y 	= geom.y
	node.SE 	= new(Node)
	newNode(node.SE,node,newGeom,3)

	//Recursively keep refining
	if refine(node.NE) {
		refineNode(node.NE)
	}

	if refine(node.NW) {
		refineNode(node.NW)
	}

	if refine(node.SW) {
		refineNode(node.SW)
	}

	if refine(node.SE) {
		refineNode(node.SE)
	}
}


/*
 * This is a helper function to construct a new leaf node
 */

func newNode(node *Node, par *Node,geom Geom, childType int) {
	node.parent 		= par
	node.currentLevel 	= par.currentLevel+1
	node.childType 		= childType
	node.isLeaf 		= true
	node.geo 		= geom
}


/*
 * This method on the segment determines if the segment and node intersect
 * Note it does not return where they intersect just if they intersect
 */

func intersect(node *Node,seg Segment) bool {
	g	:= node.geo
	minX	:= seg.x0
	maxX	:= seg.x1
	dx	:= maxX-minX

	minY	:= seg.y0
	maxY	:= seg.y1

	if maxX > (g.x + g.w) {
		maxX = g.x + g.w
	}
	if minX < g.x {
		minX = g.x
	}
	if minX >= maxX {
		return false
	} 

	if dx > 0 {
		minY = 	seg.slope * minX + seg.intercept
		maxY = 	seg.slope * maxX + seg.intercept
		if minY > maxY {
			temp	:= maxY
			maxY 	= minY
			minY 	= temp
		}
	}

	if maxY > (g.y + g.h) {
		maxY = g.y + g.h
	}
	if minY < g.y {
		minY = g.y
	}
	if minY >= maxY {
		return false
	}

	return true
}


/*
 * This function determines whether or not to refine a node
 * 
 * Refines if the level is than the maximum and 
 * the node and segment intersect
 */

func refine(node *Node) bool {
	return node.currentLevel<maxLevel && intersect(node,*segment)
}


/*
 * This functions coarsens the note by setting its children to nil
 * The children are garbage collected
 * The node is now tagged as a leaf node
 */

func coarsenNode(node *Node) {//i don't think this is how to delete?
	node.NE 	= nil
	node.NW 	= nil
	node.SW 	= nil
	node.SE 	= nil
	node.isLeaf 	= true
	node.toCoarsen 	= false
}


/*
 * This functions determines whether or to coarsen a node
 * Coarsens if the line segment does not intersect the node
 * Only called on nodes that are not leaves
 */

func coarsen(node *Node) bool {
	return !(intersect(node,*segment))
}


/*
 * This function calls the other functions that fill the neighbor slice
 * for each node
 */

func getNeighbors(node *Node, neighbors *[][]*Node) {

	n 	:= *neighbors
	if node.parent == nil {//parents have no neighbors
		return  
	}

	//Get the neighbor slice of the parent node
	parentsNeighbors := *getParentsNeighbors(node.parent)
		
	//How the neighbor slice is filled depends on its child type
	switch node.childType {
	
	case 0: //NE Child
		getNeighborsSibs(node.parent.NW,&n[3],0,3)
		getNeighborsSibs(node.parent.SE,&n[1],0,1)
		processNeighbors(parentsNeighbors[2],&n[2],node.currentLevel-1,1)
		processNeighbors(parentsNeighbors[0],&n[0],node.currentLevel-1,3)
	
	case 1: //NW Child	
		getNeighborsSibs(node.parent.NE,&n[2],1,2)
		getNeighborsSibs(node.parent.SW,&n[1],0,1)
		processNeighbors(parentsNeighbors[3],&n[3],node.currentLevel-1,0)
		processNeighbors(parentsNeighbors[0],&n[0],node.currentLevel-1,2)
	
	case 2: //SW Child
		getNeighborsSibs(node.parent.SE,&n[2],1,2)
		getNeighborsSibs(node.parent.NW,&n[0],2,3)
		processNeighbors(parentsNeighbors[3],&n[3],node.currentLevel-1,3)
		processNeighbors(parentsNeighbors[1],&n[1],node.currentLevel-1,1)
	
	case 3: //SE Child
		getNeighborsSibs(node.parent.SW,&n[3],0,3)
		getNeighborsSibs(node.parent.NE,&n[0],2,3)
		processNeighbors(parentsNeighbors[2],&n[2],node.currentLevel-1,2)
		processNeighbors(parentsNeighbors[1],&n[1],node.currentLevel-1,0)
	}
	
	*neighbors = n //update the pointer to point to the filled slice
}


/*
 * This function return a pointer to the slice of the nodes neighbors
 */
func getParentsNeighbors(node *Node) *[][]*Node {

	var s []*Node
	parentsNeighbors 		:= make([][]*Node,4)
	
	for i := 0; i < 4; i++ {
		parentsNeighbors[i]	= s
	}

	getNeighbors(node,&parentsNeighbors)
	return &parentsNeighbors
}


/*
 * This function returns the neighbors that are either the sibling of the node
 * or the children of a sibling since all neighbors are leaf nodes
 * 
 * Based on the child type of the siblings the function recurses down the tree
 * until it reaches the correct neighbor leaves and adds them to the slice
 */

func getNeighborsSibs(node *Node,list *[]*Node,sibOneType,sibTwoType int) {
	if(node.isLeaf){
		l	:= *list
		*list 	= append(l,node)
	} else {
		if sibOneType == 0 || sibTwoType == 0 {
			getNeighborsSibs(node.NE,list,sibOneType,sibTwoType)
		}
		
		if sibOneType == 1 || sibTwoType == 1 {
			getNeighborsSibs(node.NW,list,sibOneType,sibTwoType)
		}

		if sibOneType == 2 || sibTwoType == 2 {
			getNeighborsSibs(node.SW,list,sibOneType,sibTwoType)
		}

		if sibOneType == 3 || sibTwoType == 3 {
			getNeighborsSibs(node.SE,list,sibOneType,sibTwoType)
		}		
	}
}


/* 
 * This function determines whether or node is actually a neighbor
 * Based on the processing of the neighbors of the parent
 */
func (node *Node) actualNeighbors(level,nodeType int) bool {
	if node.currentLevel == level {
		return node.childType == nodeType
	} else {
		return node.parent.actualNeighbors(level,nodeType)
	}
}


/*
 * This function processes the neighbor slice of the parent node to determine
 * which of those nodes are also neighbors 
 */

func processNeighbors(pNeighbors []*Node,list *[]*Node,pLevel,nodeType int) {
	l:=*list
	if len(pNeighbors) == 0 { //parent does not have neighbors
		return
	} else if pNeighbors[0].currentLevel <= pLevel {
		//parent neighbor higher in tree than node
		*list 	= append(l,pNeighbors[0])
	} else {
		for i := range pNeighbors {
			if (pNeighbors[i].currentLevel == (pLevel+1)) &&
			(pNeighbors[i].childType == nodeType) {
				*list =	append(l,pNeighbors[i])
			} else if
			pNeighbors[i].actualNeighbors(pLevel+1,nodeType) {
				*list =	append(l,pNeighbors[i])
			}
		}
	}
}


/************************** Testing and Updating ************************/


/* 
 * This function is passed to the work queue for updating the quad tree
 */

func update(node *Node){ 
	if refine(node) {
		refineNode(node)
	} else if node.toCoarsen {
		coarsenNode(node) 
	} 
}
/*
 * This struct is passed to the work queue during concurrent update
 * It contains a node and a function that takes a node as a parameter
 */

type Work struct {
	node		*Node
	function 	func(*Node)
}


/*
 * This functions adds work to the work queue chan by traversing the tree
 */

func addWork(queue chan *Work, node *Node, f func(*Node)) {
	if node.isLeaf { 
		//all leaf nodes create work to check if they need refinement
		w 		:= new(Work)
		w.node		= node
		w.function	= f 
		queue<-w //add work
	} else if coarsen(node) { 
		//if we should coarsen the node, we send the
		//	actual coarsening work the queue
		node.toCoarsen 	= true
		w		:= new(Work)
		w.node		= node
		w.function	= f
		queue<-w //add work
	} else { //recursively add work
		addWork(queue,node.NE,f)
		addWork(queue,node.NW,f)
		addWork(queue,node.SW,f)
		addWork(queue,node.SE,f)
	}
}


/*
 * This function handles the concurrent work
 * It calls functions that
 *	create the work queue
 *	add work to queue 
 *	make sure that the work finishes
 */

func WorkParallel(node *Node, f func(*Node),done chan int) {
	queue 	:= make(chan *Work)
	ncpu	:= runtime.NumCPU() 	//number of cpus determines the number of go routines
	runtime.GOMAXPROCS(ncpu)

	spawnWorkers(queue,ncpu)	//set up the queue for the go routines
	addWork(queue,node,f)		//add work to the queue
	finished(queue,ncpu)		//all done
	done<-1				//all finished updating
}


/*
 * adds works to the queue based on the max number of goroutines 
 */

func spawnWorkers(queue chan *Work,ncpu int) {
	for i := 0; i < ncpu; i++ {
		go Worker(i,queue)
	}
}


/*
 * Execute the item on the worker queue
 */
func Worker(id int,queue chan *Work) {
	var work *Work
	for {
		work = <-queue 			//grab work from queue
		if work == nil { 		//all work done
			break
		}
		work.function(work.node)	//execute work to update the tree
	}
}


/*
 * This function passes nil to the work queue to signal that all processes are
 * complete
 */
func finished(queue chan *Work,ncpu int) {
	for i := 0; i < ncpu; i++ {
		queue<-nil
	}
}

/*
 * This struct is passed to the work queue when we find the neighbors
 * 
 * It contains a node, a function (which will be our neighbor function)
 * and the slice of neighbrs to be filled, which belong to the node
 */

type neighborWork struct {
	node 		*Node
	function 	func(*Node,*[][]*Node)
	neighbors 	[][]*Node
}


/* 
 * The functions below set up the concurrent work of neighbor finding 
 * 
 * They are very similar to the functions for concurrent update, so they do not
 * have separate documentation
 */

func addNeighborWork(queue chan *neighborWork, node *Node) {
	if node.isLeaf {
		//create work which takes a node, function, and
		//slice of slices to hold neighbors
		w		:= new(neighborWork)
		w.node		= node
		w.function	= getNeighbors

		var s []*Node
		n		:= make([][]*Node,4)
		for i := 0; i < 4; i++ {
			n[i]=s
		}
		w.neighbors 	= n
		queue<-w
	} else {
		//creat work as above and then recurse on children
		w		:= new(neighborWork)
		w.node		= node
		w.function	= getNeighbors

		var s []*Node
		n		:= make([][]*Node,4)
		for i := 0; i < 4; i++ {
			n[i]	= s
		}
		w.neighbors 	= n
		queue<-w

		addNeighborWork(queue,node.NE)
		addNeighborWork(queue,node.NW)
		addNeighborWork(queue,node.SW)
		addNeighborWork(queue,node.SE)
	}
}

func neighborWorkParallel(node *Node,done chan int) {
	queue	:= make(chan *neighborWork)

	//number of cpus determines the number of go routines 
	ncpu	:= runtime.NumCPU() 
	runtime.GOMAXPROCS(ncpu)

	spawnNeighborWorkers(queue,ncpu)	//set up ncpu number of go routines
	addNeighborWork(queue,node) 		//add work to the queue
	finishedNeighbor(queue,ncpu) 		//finished adding work
	
	//all done finding neighbors, set info back to main
	done<-1 				
}

func spawnNeighborWorkers(queue chan *neighborWork,ncpu int) {
	for i := 0; i < ncpu; i++ {
		go neighborWorker(i,queue)
	}
}

func neighborWorker(id int,queue chan *neighborWork) {
	var work *neighborWork
	for {
		work 			= <-queue 		//grad work from queue
		if work==nil {
			break 					//all done so quit
		}
		work.function(work.node,&work.neighbors) 	//get neighbors
		work.node.neighbors 	= work.neighbors 	//assign neighbors to node
	}
}

func finishedNeighbor(queue chan *neighborWork,ncpu int) {
	for i:=0;i<ncpu;i++ {
		queue<-nil
	}
}

/*
 * This struct is passed to the work queue during the timing phase
 * It simulates very time consuming work that would occur at each leaf
 * This helps simulate concurrent work
 */

 type dumbyWork struct {
	node *Node
	iter int
	function func(*Node,int)
}


/*
 * The functions below set up concurrent dumby work
 */

func dumby(node *Node,iter int) {
	for i:=0;i<iter;i++ {
		//dumby work
	}
}


/*
 * This functions executes the fake work serially 
 * Used to measure speedup
 */

func dumbyWorkSerial(node *Node,iter int) {
	if node.isLeaf {
		dumby(node,iter)
	} else {
		dumbyWorkSerial(node.NE,iter)
		dumbyWorkSerial(node.NW,iter)
		dumbyWorkSerial(node.SW,iter)
		dumbyWorkSerial(node.SE,iter)
	}
}

func addDumbyWork(queue chan *dumbyWork, node *Node, iter int,f func(*Node,int)) {
	if node.isLeaf { 	//all leaf nodes create dumby work
		w		:= new(dumbyWork)
		w.node		= node
		w.iter		=iter
		w.function	= f 
		queue<-w
	} else { 		//recursively add dumby work
		addDumbyWork(queue,node.NE,iter,f)
		addDumbyWork(queue,node.NW,iter,f)
		addDumbyWork(queue,node.SW,iter,f)
		addDumbyWork(queue,node.SE,iter,f)
	}
}

func dumbyWorkParallel(node *Node, iter int,f func(*Node,int),done chan int,cores int) {
	queue	:= make(chan *dumbyWork)
	ncpu 	:= cores 

	runtime.GOMAXPROCS(ncpu)
	spawnDumbyWorkers(queue,ncpu)		//set up the queue for the go routines
	addDumbyWork(queue,node,iter,f)		//add work to the queue
	finishedDumby(queue,ncpu)		//all done
	done<-1					//all finished updating
}

func spawnDumbyWorkers(queue chan *dumbyWork,ncpu int) {
	for i := 0; i < ncpu; i++ {
		go dumbyWorker(i,queue)
	}
}

func dumbyWorker(id int,queue chan *dumbyWork) {
	var work *dumbyWork
	for {
		work = <-queue 				//grad work from queue
		if work == nil { 				//all work done
			break
		}
		work.function(work.node,work.iter)	//execute dumby work
	}
}

func finishedDumby(queue chan *dumbyWork,ncpu int) {
	for i := 0; i < ncpu; i++ {
		queue<-nil
	}
}

/*
 * This function conducts the dumby work serially on trees of different depths
 * Used for time testing
 * The results are outputted to a csv file
 */

func depthTestSerial(level, dumbyIter,maxIter int, filename string) {
	file,err	:= os.Create(filename);
	check(err)
	l		:=[]byte("leaves,max depth,time,traversal,nodes\n")
	file.Write(l)
	

	for i := 4; i < level; i = i+2 {
		//construct tree
		maxLevel 	= i;
		root 		= Construct(nil,Geom{-4.0,-4.0,4.0,4.0},1,-1,i)
		
		//update the tree
		/*done		:=make(chan int)
		go WorkParallel(root,update,done)
		(<-done)

		//check area as sanity that update correct
		if(areaTest(root)!=(root.geo.w * root.geo.h)){
			fmt.Println("Actual Area: ",areaTest(root))
			break;
		}*/

		startTime	:=time.Now()
		for k := 0; k < maxIter; k++ {
			dumbyWorkSerial(root,0)
		}
		elapsedTime	:= time.Since(startTime)/time.Microsecond
		tTime		:=float32(elapsedTime)/float32(maxIter)

		//time the dumby work
		startTime	=time.Now()
		for j := 0; j < maxIter; j++{
			dumbyWorkSerial(root,dumbyIter)
		}
		elapsedTime	= time.Since(startTime)/time.Microsecond
		time		:= float32(elapsedTime)/float32(maxIter)
		
		//write the data to the file
		line	:=fmt.Sprintf("%d,%d,%f,%f,%d\n",countLeaves(root),i,time,
		tTime,CountNodes(root))
		l	=[]byte(line)
		file.Write(l)
	
	}	
}


/*
 * This function conducts the dumby work concurrently  
 * on trees of different depths
 * 
 * Used for time testing
 * The results are outputted to a csv file
 */

func depthTest(level,dumbyIter,maxIter int,filename string,cores int){
	file,err:= os.Create(filename);
	check(err)
	l:=[]byte("leaves,max depth,time,traversal,nodes\n")
	file.Write(l)


	for i := 4; i < level;i = i+2 {
		//construct the tree
		maxLevel = i;
		root = Construct(nil,Geom{-4.0,-4.0,4.0,4.0},1,-1,i)
		
		//update the tree
		done:=make(chan int)
		/*go WorkParallel(root,update,done)
		(<-done)

		//sanity area test for correct update
		if(areaTest(root) != (root.geo.w * root.geo.h)){
			fmt.Println("Actual Area: ",areaTest(root))
			break;
		}*/
		
		//time traversal
		
		startTime	:=time.Now()
		for k := 0; k < maxIter; k++ {
			go dumbyWorkParallel(root,0,dumby,done,cores)
			(<-done)
		}
		elapsedTime 	:= time.Since(startTime)/time.Microsecond
		tTime		:= float32(elapsedTime)/float32(maxIter)
		fmt.Println("traversal time: ", tTime)
		//timing the dumby work
		startTime	=time.Now()
		for j := 0; j< maxIter; j++{
			go dumbyWorkParallel(root,dumbyIter,dumby,done,cores)
			(<-done)
		}
		elapsedTime	= time.Since(startTime)/time.Microsecond
		time		:= float32(elapsedTime)/float32(maxIter)
		
		//output the data to the file
		line	:=fmt.Sprintf("%d,%d,%f,%f,%d\n",countLeaves(root),i,time,
		tTime,CountNodes(root))
		l	=[]byte(line)
		file.Write(l)
		
	}	
}


/************************** DEBUGGING ********************************/


/*
 * Prints the geometry of the neighbors of the node to stdout
 * Uses the helper method printNeighbors to prin the slice
 */

func printNodeNeighbors(node *Node) {
	if node.isLeaf {
		fmt.Println("Neighbors of node: ",node.geo)
		printNeighbors(node.neighbors)
	} else {
		fmt.Println("Neighbors of node: ",node.geo)
		printNeighbors(node.neighbors)
		printNodeNeighbors(node.NE)
		printNodeNeighbors(node.NW)
		printNodeNeighbors(node.SW)
		printNodeNeighbors(node.SE)
	}
}


/* 
 * Helper method to prince the geometry of the neighbors
 */

func printNeighbors(neighbors [][]*Node) {
	for i := 0; i < 4; i++ {
		switch i {
		case 0:
			fmt.Println("N neighbors")
		case 1:
			fmt.Println("S neighbors")
		case 2:
			fmt.Println("E neighbors")
		case 3:
			fmt.Println("W neighbors")
		}
		for j := 0; j < len(neighbors[i]); j++ {
			fmt.Println(neighbors[i][j].geo)
		}
	}
	fmt.Println()
}


/*
 * This function counts the total nodes in the tree
 */

func CountNodes(node *Node) int {
	if node.isLeaf {
		return 1
	}else {
		return 1 +
		CountNodes(node.NE) + CountNodes(node.NW) +
		CountNodes(node.SW) + CountNodes(node.SE)
	}
}


/*
 * This function returns the number of leaves in the tree
 */

func countLeaves(node *Node) int {
	if node.isLeaf {
		return 1;
	} else {
		return countLeaves(node.NE)+countLeaves(node.NW)+
			countLeaves(node.SW)+countLeaves(node.SE)
	}
}


/*
 * Checks if there was an error
 * Used in file i/o
 */

func check(e error) {
	if e != nil {
		fmt.Println("oops")
	    panic(e)
    }
}


/*
 * This function calculates the spatial area of the quad tree
 */

func areaTest(node *Node) float32 {
	if node.isLeaf {
		return node.geo.w * node.geo.h
	} else {
		return areaTest(node.NE)+areaTest(node.NW)+
			areaTest(node.SW)+areaTest(node.SE)
	}
}


func main() {

	//initialize segment
	originalB := float32(-7.0)
	segment =  &Segment {x0:-4.0,x1: 0.0, slope:-1.0, intercept:originalB}
	segment.y0 = segment.x0*segment.slope + segment.intercept
	segment.y1 = segment.x1*segment.slope + segment.intercept
	segment.Translate(float32(1.5),float32(1.5))//center segment	
	
	//Initialize command line arguments
	functionPtr := flag.Int("case",0,"enter number for test")
	filenamePtr := flag.String("filename","test","file for output")
	depthPtr := flag.Int("depth",0,"depth of tree")
	maxIterPtr := flag.Int("maxIter",1,"number of iterations of test")
	dumbyIterPtr := flag.Int("dumbyIter",1,"iterations of dumby work")
	numCoresPtr := flag.Int("numCores",1,"cores used")

	flag.Parse() //parse command line arguments

	//time dummy work
	node:= Construct(nil,Geom{-4.0,-4.0,4.0,4.0},1,-1,2)
	startTime	:=time.Now()
	w		:= (*dumbyIterPtr)
	for i:=0;i<(*maxIterPtr);i++ {
		dumby(node,w)
	}
	elapsedTime:=float32(time.Since(startTime)/time.Microsecond)/float32(*maxIterPtr)
	fmt.Println("Dummy Work: ",elapsedTime, "Iterations: ",w)

	switch *functionPtr {
	case 0: //concurrent depth test
		fmt.Println("Concurrent test")
		depthTest(*depthPtr,*dumbyIterPtr,*maxIterPtr,*filenamePtr,*numCoresPtr)
	case 1: //sequential depth test
		fmt.Println("Serial test")
		depthTestSerial(*depthPtr,*dumbyIterPtr,*maxIterPtr,*filenamePtr)
	case 2: //time work
		node		:= Construct(nil,Geom{-4.0,-4.0,4.0,4.0},1,-1,2)
		startTime	:= time.Now()
		w		:= (*dumbyIterPtr)
		for i:=0;i<100;i++ {
                       dumby(node,w)							                
	       	}	
	       	elapsedTime:=float32(time.Since(startTime)/time.Microsecond)/float32(100.0)
		fmt.Println("Dummy Work: ",elapsedTime,"Iterations: ",w)
	}	
}

