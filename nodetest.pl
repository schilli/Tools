#!/usr/bin/perl -w

my $VERSION = 0.1;

# ============================================================================ #

# ============= #
# USAGE MESSAGE #
# ============= #

my $USAGE = "nodetest.pl
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
    -o, --onlyempty         Use only empty nodes
    -l, --loglevel 2        Set loglevel (1 = only warnings, 2 = everything)


node names are specified with either alias or hostname indices:
      7 for alias    jubio07   (hostname: iff560c43)
    c49 for hostname iff560c49 (alias:    jubio13  )
";

# ============================================================================ #

# =============== #
# INCLUDE MODULES #
# =============== #

use warnings;
use strict;
use Getopt::Long;                  # command line parsing
use threads;                       # thread module
use threads::shared;               # enable shared memory for threads
use Term::ANSIColor qw( colored ); # for colored terminal output 
use List::Util qw(sum);            # to sum a list
use Time::HiRes qw( time );        # timing function

# ============================================================================ #

# ========== #
# PARAMETERS #
# ========== #

# Offset of jubio hostname and alias numbering
#   i. e. jubio01 => iff560c37
my $jubioOffset = 36;

# default hostnames
my @defaultHostnames = (1..32);
my $n;
foreach $n (@defaultHostnames) { $n = &hostname_from_index($n); }

# Loglevels
#my $warnings   = 1;
my $everything = 2;

# ============================================================================ #

# ==== #
# MAIN #
# ==== # 

# Global variables
my @excludeNodes = ();    # nodes not to query
my @includeNodes = ();    # nodes to query
my $qthreads     = 0;     # threads to use for node load query
my $pthreads     = 0;     # threads to use for node load parsing
my $loglevel     = 1;

my %rawInfo = ();                     # raw information about jubio node configurations and load
my %nodeUp;    share(%nodeUp);        # nodes are up (1) or down (0)
my %numCores;  share(%numCores);      # number of      cores of each jubio node
my %freeCores; share(%freeCores);     # number of free cores of each jubio node
my %memory;    share(%memory);        #               memory of each jubio node
my %freeMem;   share(%freeMem);       #          free memory of each jubio node
my %userProcs; share(%userProcs);     # running processes (or threads) per user on each jubio node (hash of hashes)

my $starttime = time();

&parse_command_line();

if ($loglevel >= $everything) {
    print "exclude:  @excludeNodes\n";
    print "include:  @includeNodes\n";
    print "qthreads: $qthreads\n";
}

&gather_jubio_info(); 

&parse_info();

&report_load();

my $endtime = time();

printf "Execution time: %.2f seconds.\n", $endtime - $starttime;

#foreach $key (keys(%rawInfo)) {print $key, " - ", @{ $rawInfo{$key} }, "\n";}
#foreach $key (keys(%rawInfo)) {print $key, "\n";}

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
        'p|pthreads=i'    => \$pthreads,
        'l|loglevel=i'    => \$loglevel);

    # print help and exit if help option was given
    if ($printhelp) {
        print $USAGE;
        exit 0;
    }

    # parse included nodes
    if (@includeNodes > 0) {
        # transform hostname or alias index into appropriate hostname
        foreach my $node (@includeNodes) { $node = &hostname_from_index($node); }
    } else {
        # if list is empty, put all standard nodes here
        @includeNodes = @defaultHostnames;
    }

    # parse excluded nodes
    if (@excludeNodes > 0) {
        # transform hostname or alias index into appropriate hostname
        foreach my $node (@excludeNodes) { $node = &hostname_from_index($node); }

        # remove every excluded node from the include array
        foreach my $node (@excludeNodes) {
            @includeNodes = grep(!/$node/, @includeNodes);
        }
    }
}

# ============================================================================ #

# convert hostname to node alias (e.g. iff560c49 -> jubio13)
sub hostname_to_alias {
    $_ = $_[0];

    if (/iff560c(\d{2})/) {
        return sprintf("jubio%02d", ($1 - $jubioOffset));

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

    my @threads;

    while ($node < @includeNodes) {
        $thread = 0;

        # create threads
        while ($thread < $qthreads and $node < @includeNodes) {
            $threads[$thread] = threads->create({'context' => 'list'}, 'gather_node_info', $includeNodes[$node]);
            $thread++;
            $node++;
        }

        # join with threads
        my $threadsCreated = $thread;
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

    my @info;

    # retrieve load information
    @info = `ssh $node "hostname; cat /proc/cpuinfo; free -m; ps aux" 2>&1`;

    #my @info = `ls -lh`;
    chomp(@info);

    # return array so we can access it in the enclosing scope
    return @info;
}

# ============================================================================ #

# parse node load information of all alive jubio nodes
sub parse_info {
    my $node = 0;
    my $thread = 0;

    my $idleCores = &idle_cores();
    if ($pthreads <= 0) {$pthreads = $idleCores;}         # if pthreads is not user specified
    if ($pthreads <= 0) {$pthreads = 1;}                  # if there are no idle cores
    
    my @threads;

    while ($node < @includeNodes) {
        $thread = 0;

        # create threads
        while ($thread < $pthreads and $node < @includeNodes) {
            $threads[$thread] = threads->create('parse_node_info', $includeNodes[$node]);
            $thread++;
            $node++;
        }

        # join with threads
        my $threadsCreated = $thread;
        for ($thread = 0; $thread < $threadsCreated; $thread++) {
            my $hostname =  $includeNodes[$node - $threadsCreated + $thread];
            $threads[$thread]->join();
        } 
    } 
}

# ============================================================================ #

sub parse_node_info {
    my ($node) = @_;

    my @info = @{ $rawInfo{$node} };

    $numCores{$node}  = 0;
    $freeCores{$node} = 0;
    $memory{$node}    = 0;
    $freeMem{$node}   = 0;

    my %procs; share(%procs);

    # if node is up
    if (not $info[0] =~ /No route to host/) {

        $nodeUp{$node} = 1;

        # loop through raw info
        foreach (@info) {

            # count processors
            if (/^processor\s+:\s+\d+$/) {
                $numCores{$node} += 1;

            # accumulate load and memory utilization
            } elsif (/^(\w+)\s+\d+\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
                $freeCores{$node} += $2;

                if (exists $procs{$1}) {
                    $procs{$1} += $2 / 100.0;
                } else {
                    $procs{$1} = $2 / 100.0;
                }

            # accumulate load and memory utilization (if ps reports load as decimal)
            } elsif (/^(\w+)\s+\d+\s+(\d+)\s+(\d+\.\d+)/) {
                 $freeCores{$node} += $2;

                if (exists $procs{$1}) {
                    $procs{$1} += $2 / 100.0;
                } else {
                    $procs{$1} = $2 / 100.0;
                } 
            }

            if (/^Mem:\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+/) {
                $memory{$node} = int($1/1024 + 0.99);
            }

            # determine free memory
            if (/^-\/\+ buffers\/cache:\s+\d+\s+(\d+)$/) {
                $freeMem{$node} = $1/1024;
            }

        }

        #$freeCores{$node} = $numCores{$node} - int($freeCores{$node} / 100.0 + 0.99);
        $freeCores{$node} = $numCores{$node} - ($freeCores{$node} / 100.0);

        $userProcs{$node} = \%procs;

    # if node is down
    } else {
        $nodeUp{$node} = 0;
    }
}

# ============================================================================ #

# estimate the number of idle cores on this machine
sub idle_cores {

    my $load = 0;
    foreach (`ps aux`) {
        if (/^\w+\s+\d+\s+(\d+\.\d+)\s+/) {
            $load += $1;
        }
    }

    my $numCores = 0;
    foreach (`cat /proc/cpuinfo`) {
        if (/^processor\s+:\s+\d+$/) {$numCores++;}
    }

    my $freeCores = $numCores - int($load/100 + 0.99);

    return $freeCores;
}

# ============================================================================ #

# report load on jubio nodes
sub report_load {

    printf "Resources reported as (free/total)\n";
    printf "%7s %9s    %5s   %6s    %s\n", "alias", "hostname", "cores", "memory", "load per user";

    my $emptyNodes     = 0;
    my $emptyNodeCores = 0;

    foreach my $hostname (sort(keys(%nodeUp))) {

        my $alias = hostname_to_alias($hostname);

        if ($nodeUp{$hostname}) {

            my $colorScheme = 'white on_black';

            if ($numCores{$hostname} - $freeCores{$hostname} < 0.1) { 
                $colorScheme = 'green on_black';
                $emptyNodes++;
                $emptyNodeCores += $freeCores{$hostname};
            }
            elsif (                    $freeCores{$hostname} < 1.0) { $colorScheme = 'red on_black'; }
            else { $colorScheme = 'yellow on_black'; }

            my $printText = sprintf "%7s %9s  %4.1f/%2d  %4.1f/%2d   ", $alias, $hostname,
                                                  $freeCores{$hostname}, $numCores{$hostname},
                                                  $freeMem{$hostname}, $memory{$hostname}; 

            print colored($printText, $colorScheme);

            my %procsPerUser = %{$userProcs{$hostname}};
            foreach my $user (keys(%procsPerUser)) {
                my $load = $procsPerUser{$user};
                if ($load > 0.1) {
                    printf " %s[%.1f]", $user, $load;
                }
            }

            print "\n";

        } else {
            my $printText = sprintf "%7s %9s     down            ", $alias, $hostname;
            print colored($printText, 'black on_red');
            print "\n";
        }
    }

    my $totalCores     = sum (values(%numCores));
    my $totalFreeCores = sum (values(%freeCores));
    my $usedCores = $totalCores - $totalFreeCores;
    my $load = $usedCores / $totalCores * 100;

    printf "Total JUBIO load: %d%% (%d of %d cores in use, %d cores free, %d empty nodes with %d free cores)\n", $load, $usedCores, $totalCores, $totalFreeCores, $emptyNodes, int($emptyNodeCores+0.5);
}

# ============================================================================ #

sub hello {
    my $ID = threads->tid();
    print "Hello World from thread $ID!\n";
    return $ID;
}
