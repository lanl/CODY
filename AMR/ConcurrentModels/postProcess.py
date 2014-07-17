#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import csv as csv
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import cPickle as pickle


data={}


# Reads in the data from the file into a dictionary
# Data keeps track of: 
	# cores (0 for serial)
	# time
	# leaves
	# traversal time
	# nodes
	# key is the depth
def readFile(filename,i,dNum):
	f=open(filename,'rt')
	reader=csv.reader(f)
	next(reader,None) #skip header
	for row in reader:
		leaves=float(row[0])
		depth=float(row[1])
		time=float(row[2])/1000000.0
		tTime=float(row[3])/1000000.0
		nodes=float(row[4])
		if(i==0):
			data[depth]=[[i],[time],leaves,[tTime],nodes]
		else:
			data[depth][0].append(i)
			data[depth][1].append(time)
			data[depth][3].append(tTime)
			#could add panic if leaves don't match
	
#Output figures for data mapped with cores on x-axis
def outputFigure(dataDic,filename,xRange,yRange,title,xLabel,yLabel,loc,
		maxCores,line):

	fig = plt.figure()
	fig.suptitle(title)
	fig.set_tight_layout(True)
	local = loc
	for key in dataDic.keys():
		local=local+1
		pt = fig.add_subplot(local)
		pt.set_xlabel(xLabel)
		pt.set_ylabel(yLabel)
		title='Depth: ' + str(key)+', Leaves: '+str(dataDic[key][2])
		pt.set_title(title)
		pt.plot(dataDic[key][0],dataDic[key][1],'bs')
		pt.axis([0,xRange,0,yRange])
		if line:
			pt.plot([0,maxCores+1],[0,maxCores+1])
		plt.tight_layout()	
	fig.savefig(filename)

def main():
	dirc = 2 #number of directories 
	outFiles = dirc * 3 #number of different dummy work times
	parser = argparse.ArgumentParser(description="command line args")
	parser.add_argument('-d','--directories',help='directory name',
			nargs=dirc,required=True)
	parser.add_argument('-m','--maxCores',type=int,help='maximum # of cores',
			required=True)
	parser.add_argument('-o','--output',help='output file name',
			nargs=outFiles,required=True)
	args = parser.parse_args()
	maxCores=args.maxCores
	outputFiles=args.output
	end='.csv'
	directories=args.directories
	indexOut = 0
	figureLocation = 220
	

	dircNum = 1
	#print nodeTime
	for directory in directories:
		start=directory+'depthTest'
		filename=start+end
		readFile(filename,0,dircNum) #serial code
		cores=1
		while cores<maxCores+1:
			filename = start+str(cores)+end
			readFile(filename,cores,dircNum)
			cores=cores*2
	#	print nodeTime
		speedUp = {}
		for key in data.keys():
			t=[]
			seq = data[key][1][0]
			for s in data[key][1]:
				t.append(seq/s)
			speedUp[key]=[data[key][0][1:],t[1:],data[key][2]]
		outputFigure(speedUp,outputFiles[indexOut],maxCores+1,maxCores+1,
		"Speed up",'number of cores','speed up',figureLocation,maxCores,
		True)
		indexOut = indexOut+1
		strongScale = {}
		for key in data.keys():
			t = []
			one = data[key][1][1]
			for s in data[key][1]:
				t.append(one/s)
			strongScale[key]=[data[key][0][1:],t[1:],data[key][2]]

		outputFigure(strongScale,outputFiles[indexOut],maxCores+1,
				maxCores+1,"Strong Scaling",'number of cores',
				'speed up',figureLocation,maxCores,True)
		
		indexOut=indexOut+1
		
		pickle.dump(data,open(outputFiles[indexOut],"wb"))	
		
		indexOut=indexOut+1
		
		dircNum = dircNum+1
		data.clear()
		speedUp.clear()
	
			
			
if __name__ == '__main__':
	main()
