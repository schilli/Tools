#!/usr/bin/perl -w

# ============================================================================ #

# ============= #
# USAGE MESSAGE #
# ============= #

$USAGE = "nodetest.pl

Autor:  Oliver Schillinger
E-Mail: o.schillinger\@fz-juelich.de

This program queries and reports the workload on jubio nodes

To be invoked as:
    nodetest.pl <options>

The following options are recognized:
    -h, --help              Print this option
    -e, --exclude 5 8 9     Exclude nodes
    -i, --include 2 7 9     Only include these nodes
    -n, --nthreads 10       Number of threads to use for node queries
                            Defaults to total number of queries

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
@defaultHostnames = (37..68);

# ============================================================================ #

# ==== #
# MAIN #
# ==== # 

# Global variables
@excludeNodes = (); # nodes not to query
@includeNodes = (); # nodes to query

&parse_command_line();

&hostname_to_alias("iff560c49");
&alias_to_hostname("jubio13");

#&get_hostname_from_alias("jubio01");

# ============================================================================ #

# =========== #
# SUBROUTINES #
# =========== #  

# ============================================================================ #

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

sub hostname_to_alias {
    my ($hostname) = @_;

    $_ = $hostname;

    if (/iff560c/) {
        my $aliasIndex = $' - $jubioOffset;
        return "jubio$aliasIndex";

    } else {

        die "Bad hostname: $hostname";

    }
}

# ============================================================================ #

sub alias_to_hostname {
    my ($alias) = @_;

    $_ = $alias;

    if (/jubio/) {
        my $jubioIndex = $' + $jubioOffset;
        return "iff560c$jubioIndex";

    } else {

        die "Bad alias: $alias";

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



