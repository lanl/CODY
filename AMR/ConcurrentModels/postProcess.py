#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import csv as csv
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import cPickle as pickle


data={}
nodeTime = {}

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
			if dNum==1:
				nodeTime[i][0][0].append(nodes)
				nodeTime[i][0][1].append(leaves)
			j=dNum+(2*(dNum-1))
			nodeTime[i][j][0].append(tTime)
			nodeTime[i][j][1].append(time)
			nodeTime[i][j+1].append(leaves/time)
			nodeTime[i][j+2].append(nodes/tTime)
		else:
			data[depth][0].append(i)
			data[depth][1].append(time)
			data[depth][3].append(tTime)
			if dNum==1:
				nodeTime[i][0][0].append(nodes)
				nodeTime[i][0][1].append(leaves)
			j=dNum+(2*(dNum-1))
			
			nodeTime[i][j][0].append(tTime)
			nodeTime[i][j][1].append(time)
			nodeTime[i][j+1].append(leaves/time)
			nodeTime[i][j+2].append(nodes/tTime)
			#could add panic if leaves don't match
	

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
	dirc = 1
	outFiles = dirc * 5
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
	
	nodeTime[0] = [[] for x in xrange(3*dirc+1)]
	nodeTime[0][0] = [[] for x in xrange(2)]
	d = 1
	while d<dirc+1:
		nodeTime[0][d+(2*(d-1))] = [[] for x in range(2)]
		d = d+1
	cores = 1
	#print "started cores"
	while cores<maxCores+1:
		nodeTime[cores]=[[] for x in xrange(3*dirc + 1)]
		nodeTime[cores][0] = [[] for x in xrange(2)]
		d = 1
		while d<dirc+1:
			nodeTime[cores][d+(2*(d-1))] = [[] for x in range(2)]
			d = d+1
	
		cores = cores*2

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
		overhead = {}
		for key in data.keys():
			t=[]
			over=[]
			seq = data[key][1][0]
			seqTime = data[key][3][0]
			for s in data[key][1]:
				t.append(seq/s)
			for o in data[key][3]:
				over.append(o-seqTime)
			speedUp[key]=[data[key][0][1:],t[1:],data[key][2]]
			overhead[key]=[data[key][0][1:],over[1:],data[key][2]]
			y = max(over)
		outputFigure(speedUp,outputFiles[indexOut],maxCores+1,maxCores+1,
		"Speed up",'number of cores','speed up',figureLocation,maxCores,
		True)
		indexOut = indexOut+1
		outputFigure(overhead,outputFiles[indexOut],maxCores+1,y,"Overhead",
		'number of cores', 'time',figureLocation,maxCores,False)
		indexOut=indexOut+1
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
		pdf_page = PdfPages(outputFiles[indexOut])
	
		for key in nodeTime.keys():
			timeFig = plt.figure()
			timeFig.suptitle(('Cores: ' + str(key)))	
			plt.xlabel('number of nodes')
			plt.ylabel('traversal time')
			plt.plot(nodeTime[key][0][0],nodeTime[key][dircNum+(2*(dircNum-1))][0],'bs')
			maxX = max(nodeTime[key][0][0])
			maxY = max(nodeTime[key][dircNum+(2*(dircNum-1))][0])
			plt.axis([0,maxX,0,maxY])
			pdf_page.savefig(timeFig)
		
		
		for key in nodeTime.keys():
			timeFig = plt.figure()
			timeFig.suptitle(('Cores: ' + str(key)))
			plt.xlabel('number of leaves')
			plt.ylabel('work time')
			plt.plot(nodeTime[key][0][1],nodeTime[key][dircNum+(2*(dircNum-1))][1],'bs')
			maxX = max(nodeTime[key][0][1])
			maxY = max(nodeTime[key][dircNum+(2*(dircNum-1))][1])
			plt.axis([0,maxX,0,maxY])
			pdf_page.savefig(timeFig)

		for key in nodeTime.keys():
			timeFig = plt.figure()
			timeFig.suptitle(('Cores: ' + str(key)))
			
			plt.xlabel('number of leaves')
			plt.ylabel('nodes per second')
			plt.plot(nodeTime[key][0][1],nodeTime[key][dircNum+(2*(dircNum-1))+1],'bs')
			maxX = max(nodeTime[key][0][1])
			maxY = max(nodeTime[key][dircNum+(2*(dircNum-1))+1])
			plt.axis([0,maxX,0,maxY])
			pdf_page.savefig(timeFig)

		for key in nodeTime.keys():
			timeFig = plt.figure()
			timeFig.suptitle(('Cores: ' + str(key)))
			
			plt.xlabel('number of nodes')
			plt.ylabel('nodes per second')
			plt.plot(nodeTime[key][0][0],nodeTime[key][dircNum+(2*(dircNum-1))+2],'bs')
			maxX = max(nodeTime[key][0][0])
			maxY = max(nodeTime[key][dircNum+(2*(dircNum-1))+2])
			plt.axis([0,maxX,0,maxY])
			pdf_page.savefig(timeFig)
	
		pdf_page.close()
	
		indexOut=indexOut+1
		pickle.dump(data,open(outputFiles[indexOut],"wb"))	
		indexOut=indexOut+1
		dircNum = dircNum+1
		data.clear()
		speedUp.clear()
	
	pickle.dump(nodeTime,open(outputFiles[indexOut],"wb"))
			
			
if __name__ == '__main__':
	main()
