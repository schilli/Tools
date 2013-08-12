#! /usr/bin/env python
#
# Tool:   running_average.py
# Author: Oliver Schillinger
# Date:   Juli 03, 2013
#
# The function running_average() computes the running average of
# a linear dataset, averaged with one of three window functions:
#   constant
#   linear
#   gaussian
#
# When called directly use as:
# running_average <outfile> <column> <width>
#   <outfile>   file to write output ot
#   <column>    column to comput running average of
#   <width>     window width
# input is read from stdin


import math, sys
import numpy as np
from scipy.constants import pi
from running_average import *


if __name__ == "__main__":
    main()  

# ============================================================================ #


def main():

    # parse cmd args
    if len(sys.argv) > 1:
        outfile = sys.argv[1]
    else:
        outfile = "average.out"

    if len(sys.argv) > 2:
        col = int(sys.argv[2]) - 1
    else:
        col = 0

    if len(sys.argv) > 3:
        width = int(sys.argv[3])
    else:
        width = 1

    # read data from stdin
    data = sys.stdin.readlines()

    # extract data
    values = np.zeros(len(data), dtype=np.float64);
    for i, line in enumerate(data):
        data[i] = line.split()
        values[i] = data[i][col]

    # compute running average
    average = running_average(values, width, "gaussian")

    # save average in file
    of = open(outfile, 'w')
    for i, line in enumerate(data):
        for j in line[:col]:
            of.write(j)
            of.write("\t")

        of.write(str(average[i]))
        of.write("\t")

        for j in line[col+1:]:
            of.write(j)
            of.write("\t")

        of.write("\n")

    of.close()

    print "Computed running average of column {} of input data ({} datapoints)".format(col+1, len(average))


# ============================================================================ #


def const_win(width):
    """return a normalized array of halfwidth width full of ones"""
    return  np.ones(2*width+1) / (2*width+1) 


# ============================================================================ #


def linear_win(width):
    """return a normalized array that linearly increases
    in the first half from 0 to 1 and decreases
    again to 0 in the second half"""

    if (width == 0):
        return np.ones(1)

    a = np.zeros(2*width+1)

    b = np.array(range(width+2)[1:], dtype=float)
    b = b / b[-1]

    a[:width+1] = b
    a[width+1:] = b[::-1][1:]

    # normalize
    a = a / sum(a)

    return a


# ============================================================================ #


def gaussian_win(width):
    """return an array that contains a normalized gaussian
    from -2*sigma to 2*sigma"""

    if (width == 0):
        return np.ones(1)

    sigma = width / 2.0
    f = 1.0/(sigma * math.sqrt(2*pi))

    a = np.array(range(-width, width+1), dtype=float)

    for i, n in enumerate(a):
        a[i] = f*math.exp(-n**2/(2*sigma**2))

    # normalize
    a = a / sum(a) 

    return a


# ============================================================================ #


def get_window(width, win="const"):
    """returns a window normalized array
    corrsponding to the function win"""

    # select window function
    if (win == "const"):
        window = const_win(width)
    elif (win == "linear"):
        window = linear_win(width)
    elif (win == "gaussian"):
        window = gaussian_win(width)
    else:
        window = const_win(width) 

    return window


# ============================================================================ #


def running_average(data, width, win="const"):
    """compute running average of data
    Data is a list of values of which to compute the
    running average. The window function to be used 
    and its half width can be specified. The result
    is returned as a numpy array.
    Options are:  const
                  linear
                  gaussain """

    L = len(data)
    average = np.zeros(L, dtype=float)

    # loop through data
    for pos1 in range(L):

        # get upper and lower indices for window
        lower = pos1 - width
        upper = pos1 + width
        window = get_window(width, win)

        if (lower < 0):
            lower = 0
            upper = 2*pos1
            window = get_window(pos1, win)

        if (upper >= L):
            lower = pos1-(L-pos1-1)
            upper = L-1
            window = get_window(L-pos1-1, win)


        # loop through window
        for i, pos2 in enumerate(range(lower, upper+1)):
            average[pos1] += data[pos2] * window[i]

    return average
