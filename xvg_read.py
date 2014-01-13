#! /usr/bin/env python

from xvg_read import *
import sys, os
import numpy as np
import running_average as ra
import matplotlib.pyplot as plt

# ============================================================================ #

if __name__ == "__main__":
    main()   

# ============================================================================ #

def main():

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        print "Provide a filename: {} <filename>".format(sys.argv[0])
        exit(1)


    data = XVG(filename)
    data.analyze()
    data.smoothen(0.01)
    data.plot()

# ============================================================================ #

class XVG:
    """Read and process xvg datasets"""

# ==================================== #

    def __init__(self, filename = None):
        """Initialize xvg data"""

        self.filename = filename

        self.nfields = 0   # number of data fields
        self.xdata   = []  # time series or similar
        self.data    = []  # list of data series as numpy arrays

        self.mean       = []  # mean per series
        self.mean_start = []  # average of per series first x percent
        self.mean_end   = []  # average of per series first x percent
        self.sd         = []  # standard deviation per series
        self.max        = []  # maximum value per series
        self.min        = []  # minimum value per series

        self.plot_title  = None
        self.plot_xlabel = None
        self.plot_ylabel = None
        self.plot_axis   = []

        if self.filename != None:
            self.read_data()

# ==================================== #

    def read_data(self):
        """read data in xvg format"""

        if not os.path.isfile(self.filename):
            print "File {} does not exist.".format(self.filename)
            exit(1)

        F = open(self.filename, 'r')
        data = F.readlines()
        F.close() 
        
        # read header
        l = 0
        line = data[l]
        while line[0] == '#' or line[0] =='@':

            if line[0] == '#':
                # comment line
                pass

            elif line[0] == '@':
                # plot parameters
                self.parse_parameter(line) 

            l += 1
            line = data[l]


        data = data[l:]
        self.nfields = len(data[0].split())


        # allocate memory
        self.data = []
        for i in range(self.nfields):
            self.data.append(np.zeros(len(data), dtype=np.float))
        self.xdata = np.zeros(len(data), dtype=np.float)


        # read data
        for l, line in enumerate(data):
            fields = line.split()
            self.xdata[l] = float(fields[0])
            for f, number in enumerate(fields[1:]):
                self.data[f][l] = float(number)


# ==================================== #

    def extract_string(self, s):
        """extract text between invertet double commas from s"""
        string = s
        first  = string.index('"')+1
        string = string[first:]
        last   = string.index('"')
        string = string[:last] 
        return string

# ==================================== #

    def parse_parameter(self, line):
        """Parse plot parameters from xvg file"""

        fields = line.split()

        if fields[1] == 'title':
            title = ' '.join(fields[2:])
            self.plot_title = self.extract_string(title)

        elif fields[1] == 'xaxis':
            if fields[2] == 'label':
                xlabel = ' '.join(fields[3:])
                self.plot_xlabel = self.extract_string(xlabel)

        elif fields[1] == 'yaxis':
            if fields[2] == 'label':
                ylabel = ' '.join(fields[3:])
                self.plot_ylabel = self.extract_string(ylabel) 

#        elif fields[1] == 'view':
#            self.plot_axis = []
#            for s in fields[2:]:
#                self.plot_axis.append(float(s.split(',')[0]))

        elif fields[1][0] == 's' and int(fields[1][1:]) >= 0:
            self.nfields = int(fields[1][1:]) + 1


# ==================================== #

    def analyze(self, p=0.5):
        """compute all kinds of statistical properties"""

        self.mean       = []  # mean per series
        self.mean_start = []  # average of per series first p percent
        self.mean_end   = []  # average of per series first p percent
        self.sd         = []  # standard deviation per series
        self.max        = []  # maximum value per series
        self.min        = []  # minimum value per series 

        splitpoint = int(0.5*len(self.data[0]))

        for d in self.data:
            self.mean.append(np.mean(d))
            self.mean_start.append(np.mean(d[:splitpoint]))
            self.mean_end.append(np.mean(d[splitpoint:]))
            self.sd.append(np.std(d, ddof=1))
            self.max.append(max(d))
            self.min.append(min(d))


# ==================================== #

    def smoothen(self, width):
        """smoothen the data with a window function average
        The window width is given in fraction of the number of datapoints"""

        if width >= 0.0 and width <= 1.0:
            ww = int(width*len(self.data[0]))
        else:
            print "ERROR: window width must be in [0.0,1.0]"
            sys.exit(1)


        for i, d in enumerate(self.data):
            self.data[i] = ra.running_average(d, ww, win="gaussian")


# ==================================== #

    def plot(self, series=[], outname=None, show=False, title='off'):
        """Plot specified datasets"""

        plt.figure()

        if len(series) == 0:
            for d in self.data:
                plt.plot(self.xdata, d)

        else:
            for s in series:
                plt.plot(self.xdata, self.data[s])

        # annotation
        if title != 'off' and self.plot_title != '':
            plt.title(self.plot_title)

        plt.xlabel(self.plot_xlabel)
        plt.ylabel(self.plot_ylabel)

        if len(self.plot_axis) > 0:
            plt.axis(self.plot_axis)


        # print to file
        if outname != None:
            plt.savefig(outname, bbox_inches='tight')

        # show plot
        if show or (not show and outname==None):
            plt.show()
