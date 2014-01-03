#!/usr/bin/perl -w

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

# ============================================================================ #

# ==== #
# MAIN #
# ==== # 

# Global variables
@excludeNodes = (); # nodes not to query
@includeNodes = (); # nodes to query

&parse_command_line();

#print &hostname_from_index(13), "\n";
#print &hostname_from_index("c49"), "\n";
#print &hostname_to_alias("iff560c49"), "\n";
#print &alias_to_hostname("jubio13")  , "\n";

# ============================================================================ #

# =========== #
# SUBROUTINES #
# =========== #  

# ============================================================================ #

# parse command line arguments and set variables accordingly
sub parse_command_line {

    my $printhelp = 0;

    GetOptions('h|help' => \$printhelp,
        'e|exclude=s{1,}' => \@excludeNodes,
        'i|include=s{1,}' => \@includeNodes);

    # print help and exit if required
    if ($printhelp) {
        print $USAGE;
        exit 0;
    }

    # parse included nodes
        # transform hostname or alias index into appropriate hostname
        # if list is empty, put all standard nodes here
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

sub get_hostname_from_alias {
    my $alias = $_[0];

    @data = `ssh $alias "hostname; hostname -a"`;

    print $data[0];
    print $data[1];
}

# ============================================================================ #

#$workmax  = 1000;
#$nthreads = 2;
#
#for ($i = 0; $i < $nthreads; $i++) {
#    $threads[$i] = threads->create('work');
#}
#
#for ($i = 0; $i < $nthreads; $i++) {
#    $thrReturn[$i] = $threads[$i]->join();
#    print "Thread $i returned with $thrReturn[$i]\n";
#}
#
#
#sub hello {
#    $ID = threads->tid();
#    print "Hello World from thread $ID!\n";
#}
#
#sub work {
#    my $i;
#    my $j;
#    for ($i = 0; $i < $workmax; $i++) {
#        for ($j = 0; $j < $workmax; $j++) {
#        }
#    }
#    return $i * $j;
#}



