#!/usr/bin/env python

import sys
import matplotlib.pyplot as pyplot
import numpy

class Vol2birdVisualizer():

    def __init__(self,theFile,nLinesSkip):

        self.theFile = theFile
        if nLinesSkip < 0:
            print "second input argument should be positive integer number"
            sys.exit(1)
        
        self.nLinesSkipMeta = nLinesSkip - 10
        self.nLinesSkip = nLinesSkip
        self.meta = dict.fromkeys(["heig","elev","nRang","nAzim","rangeScale","azimScale","valueOffset","valueScale","missing","product"])
        self.data = []
        
        f = open(self.theFile, 'r')
        
        self.skipLines(f)
        self.readMetaData(f)
        self.readData(f)
        
        f.close()
        
        self.plot()


    def skipLines(self,f):
        
        # skip self.nLinesSkipMeta lines
        for k in xrange(0,self.nLinesSkipMeta):
            line = f.readline()

    def readMetaData(self,f):
        
        line = f.readline().split()
        varName = line[0]
        value = line[2]
        if varName == "product":
            self.meta["product"] = value
        else:
            print "Something's not right"
        
        
        nMetaLines = 9
        for iMetaLine in xrange(0,nMetaLines):
            line = f.readline()
            lineList = line.split()
            varName = lineList[0].split('->')[1]
            value = lineList[2]
            if varName in ["heig","elev","rangeScale","azimScale","valueOffset","valueScale"]:
                self.meta[varName] = float(value)
            elif varName in ["nRang","nAzim","missing"]:
                self.meta[varName] = int(value)
            else:
                print "Something's not right"
                
                
        
    def readData(self,f):
        
        nLines = self.meta["nAzim"]
        for iLine in xrange(0,nLines):
            line = f.readline()
            lineList = [int(x) for x in line.split()]
            self.data.append(lineList)
            
            
    def plot(self):
 
        iAzim = xrange(0,self.meta["nAzim"])
        iRang = xrange(0,self.meta["nRang"])

        #setup the 2D grid with Numpy
        rang, azim = numpy.meshgrid(iRang, iAzim)

        #convert list of lists 'self.data' to a numpy array for plotting
        data = numpy.array(self.data)
        
        #now just plug the data into pcolormesh, it's that easy!
        pyplot.pcolormesh(rang, azim, data)
        if len(numpy.unique(data)) == 1:
            pass
            #pyplot.colorbar([-0.5,0.5])
        else:
            pyplot.colorbar()
        pyplot.autoscale(enable=True, axis=u'both', tight=u'both')
        pyplot.title("elevation angle = " + str(self.meta["elev"]) + " degrees; product = "+ self.meta["product"])
        pyplot.xlabel("iRang")
        pyplot.ylabel("iAzim")
        pyplot.show() 
        

if __name__ == "__main__":


    def printHelp():

        print
        print "Here's how to use " + sys.argv[0]
        print
        print "Make sure some of the PRINT_* options in 'options.conf' are enabled."
        print 
        print "Run 'baltradtest.exe' and write the output to a file, using e.g"
        print "./testbaltrad.exe ../data/T_PAGZ51_C_EHDB_20141219134510.h5 2> testbaltrad.log"
        print 
        print "Open the logfile and check at which line a given block of data begins."
        print 
        print "For visualizing a block of data that starts at line 26, issue:"
        print sys.argv[0] + " testbaltrad.log 26"
        print
        print "The arguments are:"
        print "1. the log file containing the stderr output from testbaltrad.exe;"
        print "2. the line number where the data block of interest starts"
        print



    nArgs = len(sys.argv)
    
    if nArgs == 2 and sys.argv[1] in ["-h", "--help"]:

        printHelp()
        sys.exit(-1)
    elif nArgs == 3:
        pass
    else:
        printHelp()
        sys.exit(1)

    theFile = sys.argv[1]

    try:
        nLinesSkip = int(sys.argv[2]) - 1
    except ValueError:
        print "Second input argument is not an integer"
        sys.exit(1)
    except:
        print "an error occurred"
        sys.exit(1)

    vol2birdVisualizer = Vol2birdVisualizer(theFile,nLinesSkip)


