#!/usr/bin/perl -w

use threads;
use threads::shared;

# ============================================================================ #

# ========== #
# PARAMETERS #
# ========== #

# Offset of jubio hostname and alias numbering
#   i. e. jubio01 => iff560c37
#$jubioOffset = 36;

# create hashes for alias to hostname conversion

# ============================================================================ #

# ==== #
# MAIN #
# ==== # 

&get_hostname_from_alias("jubio01");

# ============================================================================ #

# =========== #
# SUBROUTINES #
# =========== #  

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



