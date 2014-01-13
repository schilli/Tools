#!/usr/bin/perl -w
# ToDo: Retrieve raw info

$VERSION = 0.1;

# ============================================================================ #

# ============= #
# USAGE MESSAGE #
# ============= #

$USAGE = "nodetest.pl
version $VERSION

Autor:  Oliver Schillinger
E-Mail: o.schillinger\@fz-juelich.de

This program queries and reports the workload on jubio nodes.
If a command is provided it will be executed on free jubio nodes.

To be invoked as:
    nodetest.pl <options>

The following options are recognized:
    -h, --help              Print this option
    -e, --exclude 5 8 9     Exclude nodes
    -i, --include 2 7 9     Only include these nodes
    -q, --qthreads 10       Number of threads to use for node queries
                            Defaults to total number of queries
    -p, --pthreads 4        Number of threads to use for node data parsing
                            Defaults to the estimated number of idle cores
    -c, --cmd <command>     Command to execute on jubio as:
                            mpiexec -n <n> <command>
    -n, --nprocs <n>        Number of processes to start
    --nohup <filename>      Filename for nohup output
                            If not specified, nohup is not used.
    -m, --mpd               Switch to start a new mpd ring on the requested hosts
                            Any ring already running will be destroyed
    -j, --interactive       Switch to specify jubio nodes to run on interactively
    -l, --loglevel 2        Set loglevel (1 = only warnings, 2 = everything)


node names are specified with either alias or hostname indices:
      7 for alias    jubio07   (hostname: iff560c43)
    c49 for hostname iff560c49 (alias:    jubio13  )
";

# ============================================================================ #

# =============== #
# INCLUDE MODULES #
# =============== #

use Getopt::Long;     # command line parsing
use threads;          # thread module
use threads::shared;  # enable shared memory for threads

# ============================================================================ #

# ========== #
# PARAMETERS #
# ========== #

# Offset of jubio hostname and alias numbering
#   i. e. jubio01 => iff560c37
$jubioOffset = 36;

# default hostnames
@defaultHostnames = (1..32);
foreach $n (@defaultHostnames) { $n = &hostname_from_index($n); }

# Loglevels
#$warnings   = 1;
$everything = 2;

# ============================================================================ #

# ==== #
# MAIN #
# ==== # 

# Global variables
@excludeNodes = ();    # nodes not to query
@includeNodes = ();    # nodes to query
$qthreads     = 0;     # threads to use for node load query
$loglevel     = 1;

%rawInfo = ();         # raw information about jubio node configurations and load
share(%numCores);      # number of      cores of each jubio node
share(%freeCores);     # number of free cores of each jubio node
share(%memory);        #               memory of each jubio node
share(%freeMem);       #          free memory of each jubio node
share(%userProcs);     # running processes (or threads) per user on each jubio node (hash of hashes)


&parse_command_line();

if ($loglevel >= $everything) {
    print "exclude:  @excludeNodes\n";
    print "include:  @includeNodes\n";
    print "qthreads: $qthreads\n";
}

&gather_jubio_info(); 

foreach $key (keys(%rawInfo)) {print $key, " - ", @{ $rawInfo{$key} }, "\n";}

# ============================================================================ #

# =========== #
# SUBROUTINES #
# =========== #  

# ============================================================================ #

# parse command line arguments and set variables accordingly
sub parse_command_line {

    my $printhelp = 0;

    GetOptions('h|help'   => \$printhelp,
        'i|include=s{1,}' => \@includeNodes,
        'e|exclude=s{1,}' => \@excludeNodes,
        'q|qthreads=i'    => \$qthreads,
        'l|loglevel=i'    => \$loglevel);

    # print help and exit if help option was given
    if ($printhelp) {
        print $USAGE;
        exit 0;
    }

    # parse included nodes
    if (@includeNodes > 0) {
        # transform hostname or alias index into appropriate hostname
        foreach $node (@includeNodes) { $node = &hostname_from_index($node); }
    } else {
        # if list is empty, put all standard nodes here
        @includeNodes = @defaultHostnames;
    }

    # parse excluded nodes
    if (@excludeNodes > 0) {
        # transform hostname or alias index into appropriate hostname
        foreach $node (@excludeNodes) { $node = &hostname_from_index($node); }

        # remove every excluded node from the include array
        foreach $node (@excludeNodes) {
            @includeNodes = grep(!/$node/, @includeNodes);
        }
    }
}

# ============================================================================ #

# convert hostname to node alias (e.g. iff560c49 -> jubio13)
sub hostname_to_alias {
    $_ = $_[0];

    if (/iff560c(\d{2})/) {
        return "jubio" . ($1 - $jubioOffset);

    } else {
        die "Bad hostname: $_";
    }
}

# ============================================================================ #

# convert node alias to hostname (e.g. jubio13 -> iff560c49)
sub alias_to_hostname {
    $_ = $_[0];

    if (/jubio(\d{2})/) {
        return "iff560c" . ($1 + $jubioOffset);

    } else {
        die "Bad alias: $_";
    }
}

# ============================================================================ #

# get a hostname from either jubio alias index /\d{2}/ or hostname index /c\d{2}/
sub hostname_from_index {

    $_ = $_[0];

    if (/(c\d{2})/) { # if it's a jubio hostname index
        return "iff560" . $1;

    } elsif (/(\d{2})/) { # if it's a jubio alias index
        return alias_to_hostname("jubio" . $1); 

    } elsif (/(\d{1})/) { # if it's a jubio alias index with only one digit
        return alias_to_hostname("jubio0" . $1);  

    } else {
        die "Bad alias or hostname index: $_";
    }

}

# ============================================================================ #

# gather jubio load information
sub gather_jubio_info {
    my $node = 0;
    my $thread = 0;
    if ($qthreads == 0) {$qthreads = @includeNodes;}

    while ($node < @includeNodes) {
        $thread = 0;

        # create threads
        while ($thread < $qthreads and $node < @includeNodes) {
            $threads[$thread] = threads->create({'context' => 'list'}, 'gather_node_info', $includeNodes[$node]);
            $thread++;
            $node++;
        }

        # join with threads
        $threadsCreated = $thread;
        for ($thread = 0; $thread < $threadsCreated; $thread++) {
            my $hostname =  $includeNodes[$node - $threadsCreated + $thread];
            @{$rawInfo{$hostname}} = $threads[$thread]->join();
        } 
    }
}
 
# ============================================================================ #

# gather information of one specific jubio node
sub gather_node_info {
    my ($node) = @_;

    # retrieve load information
#    my @info = `ssh $node "hostname; uname" 2>&1`;
    my @info = `ls -lh`;
    chomp(@info);

    # return array so we can access it in the enclosing scope
    return @info;
}

# ============================================================================ #

sub hello {
    my $ID = threads->tid();
    print "Hello World from thread $ID!\n";
    return $ID;
}
