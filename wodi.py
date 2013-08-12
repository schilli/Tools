#! /usr/bin/env python
#
# Work Distributor
# Author: Oliver Rusche
# E-Mail: o.rusche@fz-juelich.de
#
# Determine workload on the jubio cluster and
# distribute the given job as balanced as possible
# across the given number of free cores
#
# Usage:
#   wodi2.py -h | --help for usage information
#
# supply file names and such in the Settings class

from wodi2 import *
import sys, argparse, os, shutil, time

if __name__ == "__main__":
    main()  


# ============================================================================ #


def main():

    settings = Settings()

    print settings
    print ""

    # query node load
    print "Gathering information about jubio load."
    nodes = os.system(settings.nodetestpath + "> " +  settings.nodetmpfile)

    # read node load from temporary file
    nodetmp = open(settings.nodetmpfile, 'r')
    load = nodetmp.readlines()
    nodetmp.close()
    os.remove(settings.nodetmpfile)

    # list nodes and free cores
    node       = []
    freecores  = []
    totalcores = []
    for host in load:
        host = host.split()
        if (int(host[4]) > 0 and "FREE" in host):
            node.append(host[0][-2:])
            freecores.append(int(host[4]))
            totalcores.append(int(host[7]))

    # sort nodes according to free cores
    free  = sorted(zip(freecores, node), reverse=True)
    total = sorted(zip(freecores, totalcores), reverse=True)
    freecores1, node       = zip(*free)
    freecores2, totalcores = zip(*total)
    freecores1  = list(freecores1)
    freecores2  = list(freecores2)
    freecores  = freecores1
    node       = list(node)
    totalcores = list(totalcores)

    # to make sure processes don't block each other, subtract savety margin if not all cores are free
    indeedFree = list(freecores) # cores that are indeed free, ignoring safety margin
    if settings.get_safety():
        for i, free in enumerate(freecores):
            if free < totalcores[i]:
                freecores[i] = free - 1

    # check if enough nodes are available
    if sum(freecores) < settings.get_nc():
        print "Not enough cores available (only", sum(freecores), "free)" 
        sys.exit(1)

    # set up initial work distribution
    load = []
    while sum(load) < settings.get_nc():
        load.append(freecores[len(load)])
    load[-1] = settings.get_nc() - sum(load[:-1])

    # balance workload
    load = balance(load, freecores, "max")

    # report load balance
    print "The following configuration will be run on jubio:"
    print "= cores requested"
    print "- cores free"
    print "* cores in use"
    for i, l in enumerate(load):
        print "jubio{}: {:>2} cores {}{}{}".format(node[i], l, l*"=", (indeedFree[i]-l)*"-", (totalcores[i]-indeedFree[i])*"*")

    # if we must set up mpd ring ourselves,
    # check if we are on one of the nodes we want to run on
    if settings.get_mpd():

        hostnum = hostname_to_node(os.uname()[1])
        if hostnum not in node[:len(load)]:
            print "\nPlease log into one of the nodes above"
            sys.exit(0) 

        # exit running mpdring
        print "\nExiting mpd ring if any."
        os.system("mpdallexit &> /dev/null")

        # set up new hosts file
        hosts = open(settings.tmphosts, 'w')
        for i in range(len(load)):
            hosts.write("jubio"+node[i]+"\n")
        hosts.close()

        # start new mpdring
        print "Starting new mpd ring."
        #print "mpdboot -r ssh -n {} -f {}".format(len(load), settings.tmphosts)
        os.system("mpdboot -r ssh -n {} -f {}".format(len(load), settings.tmphosts))
        os.remove(settings.tmphosts)

    # construct command line for mpi execution
    exe = settings.get_program() + " "
    for arg in settings.get_arguments():
        exe += arg + " " 
    command = "mpiexec -n " + str(load[0]) + " -host " + node_to_hostname(node[0]) + " " + exe

    if len(load) > 1:
        command += ": \\\n"

    for i, l in enumerate(load[1:-1]):
        command += "        -n " + str(l) + " -host " + node_to_hostname(node[i+1]) + " " + exe + ": \\\n"
    if len(load) > 1:
        command += "        -n " + str(load[-1]) + " -host " + node_to_hostname(node[len(load)-1]) + " " + exe

    print "\nThe following command will be executed:"
    print command 

    # ask user for confirmation
    if settings.interactive:
        print "\nDo you want to proceed?"
        if raw_input("(y/n): ") != 'y':
            sys.exit(0)
        print ""


    # flush buffered output
    sys.stdout.flush()

    # run command
    if settings.get_nohup():
        if os.path.isfile("nohup.out"):
            counter=1
            while os.path.isfile("#nohup.out.{}#".format(counter)):
                counter += 1
            shutil.copyfile("nohup.out", "#nohup.out.{}#".format(counter))
            os.remove("nohup.out")
            print "Backed up file 'nohup.out' to '#nohup.out.{}#'".format(counter)
        os.system("nohup " + command + " &")

        # wait for nohup to report message
        time.sleep(1.0) 

        print "Job will run in background until it finishes or crashes."
        print "You can now log out." 

    else:
        os.system(command)




# ============================================================================ #


# translate jubio node number to hostname
def node_to_hostname(node):
    return "iff560c" + str(36 + int(node))


# ============================================================================ #


# translate hostname to jubio node number as a string
def hostname_to_node(hostname):
    if hostname == "iff560":
        return "00"
    else:
        return "{:0>2}".format(int(hostname[-2:])-36)
    

# ============================================================================ #


# Balance workload
def balance(load, freecores, mode="balanced"):
    """Balance workload for jubio
    modes:
        balanced: Make processes per host as equal as possible among nodes
        max:      Maximize the fully loaded nodes"""

    if mode == "balanced":
        # very primitive routine to achieve approximately equal balance
        while (freecores[len(load)-1] - load[-1] > 0 and load[-1] != max(load)):
            load[load.index(max(load))] -= 1
            load[-1] += 1 

    elif mode == "max":
        # balance only the last two nodes
        if len(load) > 1:
            lastTwo = load[-2] + load[-1]
            load[-2] = (lastTwo+1)/2
            load[-1] = lastTwo - load[-2]

    else:
        print "Workload balance mode unknown."
        sys.exit(1)

    return load


# ============================================================================ #


class Settings:
    """Settings and options are parsed and stored"""

    interactive = True
    nc  = 1        # number of cores to use (processes to start)
    mpd = True     # whether to start mpd ring or use existing one
    nohup = True   # whether to use nohup
    safety = True  # whether to keep a safety margin on workload to prevent job interference
    program = None # program to start
    args    = None # artguments to pass to program

    # path to nodetest script
    nodetestpath = "/usr/users/iff_th2/oschill/nodetest.sh"
    
    # some temporary files
    nodetmpfile  = "/usr/users/iff_th2/oschill/node.tmp"
    tmphosts     = "/usr/users/iff_th2/oschill/mpd.hosts.tmp" 

# ==================================== #

    def __init__(self):
        """Initialize object"""

        self.parse_cmd_arguments()

# ==================================== #

    def parse_cmd_arguments(self):
        """Parse command line arguments"""

        parser = argparse.ArgumentParser(description='Work Distributor for JUBIO cluster.')

        parser.add_argument('-n', type=int, default=1, help='Number of cores requested for the job') 
        parser.add_argument('--noi', action='store_false', help='Switch interactive mode off') 
        parser.add_argument('--mpd', default="on", choices=["on", "off"], help='Whether to start mpd ring or use existing one')
        parser.add_argument('--nohup', default="on", choices=["on", "off"], help='Whether to use nohup')
        parser.add_argument('--safety', default="on", choices=["on", "off"], help='Safety margin on node workload')
        parser.add_argument('--cmd', nargs=argparse.REMAINDER, default=[], help='Program to start on machines plus arguments')

        args = parser.parse_args()        

        # interactive mode
        self.interactive = args.noi

        # number cores
        self.nc  = args.n
        
        # mpd ring on or off
        if args.mpd == "on":
            self.mpd = True
        else:
            self.mpd = False

        # nohup on or off
        if args.nohup == "on":
            self.nohup = True
        else:
            self.nohup = False 

        # safety margin on or off
        if args.safety == "on":
            self.safety = True
        else:
            self.safety = False

        # program and arguments
        try:
            self.program = args.cmd[0]
        except IndexError:
            parser.parse_args(["-h"])
            
        self.arguments = args.cmd[1:]
        
# ==================================== #

    def get_nc(self):
        """Query number of cores to use (processes to start)"""
        return self.nc

# ==================================== #

    def get_mpd(self):
        """Query whether to start mpd ring or use existing one"""
        return self.mpd

# ==================================== #

    def get_nohup(self):
        """Query whether to use nohup"""
        return self.nohup

# ==================================== #

    def get_safety(self):
        """Query whether to use safety margin"""
        return self.safety
 
# ==================================== #

    def get_program(self):
        """Query program to start"""
        return self.program

# ==================================== #

    def get_arguments(self):
        """Query arguments to be passed to program"""
        return self.arguments

# ==================================== #

    def __str__(self):
        """Return a string of settings and options"""

        string  = "Settings:\n"
        string += "\tnc        = {}\n".format(self.nc)
        string += "\tnohup     = {}\n".format(self.nohup)
        string += "\tmpd       = {}\n".format(self.mpd)
        string += "\tsafety    = {}\n".format(self.safety)
        string += "\tprogram   = {}\n".format(self.program)
        string += "\targuments = "

        for arg in self.arguments:
            string += "{} ".format(arg)

        return string
        
# ============================================================================ #
