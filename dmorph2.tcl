# The code between opening and closing <plagiarism></plagiarism> was taken from morph.tcl from the VMD scripts library

# ============================================================================ #

# <plagiarism> 
proc morph_linear {t N} {
  return [expr {double($t) / double($N)}]
}
proc morph_cycle {t N} {
  global M_PI
  return [expr {(1.0 - cos( $M_PI * double($t) / ($N + 1.0)))/2.0}]
}
proc morph_sin2 {t N} {
  global M_PI
  return [expr {sqrt(sin( $M_PI * double($t) / double($N) / 2.0))}]
}
# </plagiarism> 

# ============================================================================ #

proc angle_diff {angle_init angle_final} {
    set angle_diff {}
    foreach ai $angle_init af $angle_final {
        set d [expr $af - $ai]
        while {$d >  180} {set d [expr $d - 360.0]}
        while {$d < -180} {set d [expr $d + 360.0]}
        lappend angle_diff $d
    }
    return $angle_diff
}

# ============================================================================ #

proc interpolate_dihedrals {molid nf f resnums phi_init psi_init phi_diff psi_diff} {
    set dphi [vecscale $f $phi_diff]
    set dpsi [vecscale $f $psi_diff]

    set phi  [vecadd $phi_init $dphi]
    set psi  [vecadd $psi_init $dpsi]

    # wrap dihedrals back in valid range
    set i 0
    foreach a $phi b $psi {
        while {$a >  180} {set a [expr $a - 360.0]}
        while {$a < -180} {set a [expr $a + 360.0]} 
        while {$b >  180} {set b [expr $b - 360.0]}
        while {$b < -180} {set b [expr $b + 360.0]}  
        set phi [lreplace $phi $i $i $a]
        set psi [lreplace $psi $i $i $b]
        incr i
    }

    # set dihedral angles
    foreach resnum $resnums a $phi b $psi {
        set r [atomselect $molid "residue $resnum" frame $nf]
        $r set phi $a
        $r set psi $b
    }
}

# ============================================================================ #
 
proc dmorph {molid N {morph_type morph_linear}} {
    # <plagiarism>
    # make sure there are only two animation frames
    if {[molinfo $molid get numframes] != 2} {
	    error "Molecule $molid must have 2 animation frames"
    }
    # workaround for the 'animate dup' bug; this will translate
    # 'top' to a number, if needed
    set molid [molinfo $molid get id]

    # Do some error checking on N
    if {$N != int($N)} {
	    error "Need an integer number for the number of frames"
    }
    if {$N <= 2} {
	    error "The number of frames must be greater than 2"
    } 

    # determine initial and final dihedral angles
    set alpha     [atomselect $molid "name CA" frame 0]
    set resnums   [$alpha get residue]
    set phi_init  [$alpha get phi]
    set psi_init  [$alpha get psi]
    $alpha frame 1
    set phi_final [$alpha get phi]
    set psi_final [$alpha get psi] 
    set phi_diff  [angle_diff $phi_init $phi_final]
    set psi_diff  [angle_diff $psi_init $psi_final]

    puts "psi0: [lindex $psi_init 0] -> [lindex $psi_final 0] (diff = [lindex $psi_diff 0])"

    # Make N-2 new frames (copied from the last frame)
    for {set i 2} {$i < $N} {incr i} {
	    animate dup frame 1 $molid
    }
    # there are now N frames
    # </plagiarism>
    
    # select atoms for superposition
    set all         [atomselect $molid "all"]
    set backbone    [atomselect $molid "backbone"]
    set backbone_t0 [atomselect $molid "backbone" frame 0]
    
#    set sumf 0
   
    # <plagiarism>
    # Do the linear interpolation in steps of 1/N so
    # f(0) = 0.0 and f(N-1) = 1.0
    for {set t 1} {$t < [expr $N - 1]} {incr t} {

#        puts -nonewline [format "\rProgress: %5.1f%%" [expr 100.0*$t/($N-2)]]
#        flush stdout

	    # Here's the call to the user-defined morph function
	    set f [$morph_type $t [expr $N - 2]]; # f grows from 0 to 1
        # </plagiarism>
        puts "$t $f"
#	    set fprev [$morph_type [expr $t - 1] [expr $N - 2]]; # f grows from 0 to 1
#        set df   [expr $f - $fprev]
#        set sumf [expr $sumf + $df]
#        set f [expr 1.0 / ($N - 1 - $t)]

        # copy coordinates from previous frame
        set currentFrame  [atomselect $molid "all" frame $t]
        set previousFrame [atomselect $molid "all" frame [expr $t - 1]]
        $currentFrame set {x y z} [$previousFrame get {x y z}]
 
        # set new dihedrals
        interpolate_dihedrals $molid $t $f $resnums $phi_init $psi_init $phi_diff $psi_diff
        $alpha frame $t
#        set psi0 [lindex [$alpha get psi] 0]
#        puts "psi0: $psi0"

        
        # superimpose on previous (first) frame
        $all frame $t
        $backbone frame $t
        $backbone_t0 frame [expr $t - 1]
        set M [measure fit $backbone $backbone_t0]
        $all move $M
    }

    puts $t

    $all frame $t
    $backbone frame $t
    $backbone_t0 frame [expr $t - 1]
    set M [measure fit $backbone $backbone_t0]
    $all move $M 

    foreach rn $resnums phi_i $phi_init psi_i $psi_init phi_f $phi_final psi_f $psi_final phi_d $phi_diff psi_d $psi_diff {
        set output [format "" $rn $phi_i $phi_f $phi_d $psi_i $psi_f $psi_d]
        puts "$rn $phi_i\t$phi_d \t$psi_d"
    }
}

