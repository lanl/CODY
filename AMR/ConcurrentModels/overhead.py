from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import csv as csv
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import cPickle as pickle


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
		title='Depth:'+ str(key)+',Leaves: '+str(dataDic[key][2])
            	pt.set_title(title)
	    	pt.plot(dataDic[key][0],dataDic[key][1],'bs')
	    	pt.axis([0,xRange,0,yRange])
	    	if line:
			pt.plot([0,maxCores+1],[0,maxCores+1])
			plt.tight_layout()
	fig.savefig(filename)

def main():
	parser = argparse.ArgumentParser(description="command line args")
	inputFiles = 4
	outputFiles = 4
	parser.add_argument('-i','--inputF',help='input file name',required=True)
	parser.add_argument('-o','--output',help='output file name',required=True)
	args = parser.parse_args()
	inputFile = args.inputF
	outputFile = args.output

	for f in range(1,inputFiles+1):
		inFile = inputFile+str(f)+".p"
		outFile = outputFile+str(f)+".png"

		data = pickle.load(open(inFile,"rb"))
		overhead = {}
		figureLocation = 220
		maxCores = 32
		for key in data.keys():
			seqTime = data[key][3][0]
			o = [x-seqTime for x in data[key][3]]
			over = [(x/y)*100 for x,y in zip(o,data[key][1])]
			overhead[key] = [data[key][0][1:],over[1:],data[key][2]]
		outputFigure(overhead,outFile,maxCores+1,m,"Overhead",'number of cores',
				'percent overhead',figureLocation,maxCores,False)
		data.clear()
		overhead.clear()

if __name__ == '__main__':
	main()
