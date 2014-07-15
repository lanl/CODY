#!/bin/bash

cd goTree
go build -o QuadTree

func=-1 #which test am I running initialized to invalid value
testName='depthTest' #data output file header
end='.csv' #data output file extension
maxDepth=11 #maximum depth of the tree
maxIter=1 #number of times to run the work queue
dumbyIter=$1 #number of of dummy work iterations (command line argument)

#sequential program
func=1
file=$testName$end
./QuadTree -case=$func -filename=$file -depth=$maxDepth -maxIter=$maxIter -dumbyIter=$dumbyIter -numCores=1

sleep 10

#concurrent program
func=0
COUNTER=1 #
maxCores=33 #one more than max cores used, cores increase by powers of two

while [ $COUNTER -lt $maxCores ]; do
	file=$testName$COUNTER$end
	./QuadTree -case=$func -filename=$file -depth=$maxDepth -maxIter=$maxIter -dumbyIter=$dumbyIter -numCores=$COUNTER
	let COUNTER=COUNTER*2
	sleep 10
done

