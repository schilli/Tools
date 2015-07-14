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

proc get_dihedrals {molID framenum} {
    set CAsel [atomselect $molID "name CA" frame $framenum]
    set nres [llength [$CAsel get residue]]

    set phi [list]
    set psi [list]

    set psisel [atomselect $molID "(residue 0 and name N CA C) or (residue 1 and name N)" frame $framenum]
    lappend psi [measure dihed [$psisel get index] frame $framenum]

    for {set i 1} {$i < [expr $nres - 1]} {incr i} {
        set phisel [atomselect $molID "(residue [expr $i - 1] and name C) or (residue $i and name N CA C)" frame $framenum] 
        set psisel [atomselect $molID "(residue $i and name N CA C) or (residue [expr $i  + 1] and name N)" frame $framenum]
        lappend phi [measure dihed [$phisel get index] frame $framenum]
        lappend psi [measure dihed [$psisel get index] frame $framenum]
    }

    set phisel [atomselect $molID "(residue [expr $nres - 2] and name C) or (residue [expr $nres - 1] and name N CA C)" frame $framenum] 
    lappend phi [measure dihed [$phisel get index] frame $framenum]

    set dihedrals [list $phi $psi]
    return $dihedrals

    # to get the values out of the function:
    # lassign [get_dihedrals top] phi psi
}

# ============================================================================ #
 

proc get_rotation_vector {molID frame1 frame2} {
    # get the vector of angles in degrees that will rotate the dihedrals of frame1 into that of frame2
    lassign [get_dihedrals $molID $frame1] phi_first psi_first
    lassign [get_dihedrals $molID $frame2] phi_final psi_final
    set phi_vec [list]
    set psi_vec [list]

    for {set i 0} {$i < [llength $phi_first]} {incr i} {
        set phi_diff [expr [lindex $phi_final $i] - [lindex $phi_first $i]]
        set psi_diff [expr [lindex $psi_final $i] - [lindex $psi_first $i]]
        if {$phi_diff >  180} {set phi_diff [expr 1 * ($phi_diff - 360.0)]}
        if {$psi_diff >  180} {set psi_diff [expr 1 * ($psi_diff - 360.0)]}
        if {$phi_diff < -180} {set phi_diff [expr 1 * ($phi_diff + 360.0)]}
        if {$psi_diff < -180} {set psi_diff [expr 1 * ($psi_diff + 360.0)]}

#        # invert
#        set phi_diff [expr -1 * $phi_diff / abs($phi_diff) * (360 - abs($phi_diff))]
#        set psi_diff [expr -1 * $psi_diff / abs($psi_diff) * (360 - abs($psi_diff))]

        lappend phi_vec $phi_diff
        lappend psi_vec $psi_diff
    }
    
    set rotation_vec [list $phi_vec $psi_vec]
    return $rotation_vec
}

# ============================================================================ #

proc rotation_matrix_arbitrary_line {point direction angle} {
    # compute the 4D rotatin matrix around an arbitrary line defined by a point the
    # line goes through and a direction vector by an angle given in degree
    # Source:
    # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    global M_PI
    set dir [vecnorm $direction]
    lassign $point a b c
    lassign $dir   u v w
    set theta [expr $M_PI * $angle / 180.0]

    set M11 [expr $u**2 + ($v**2 + $w**2) * cos($theta)]
    set M22 [expr $v**2 + ($u**2 + $w**2) * cos($theta)]
    set M33 [expr $w**2 + ($u**2 + $v**2) * cos($theta)]

    set M12 [expr $u*$v*(1-cos($theta)) - $w*sin($theta)]
    set M13 [expr $u*$w*(1-cos($theta)) + $v*sin($theta)]
    set M21 [expr $u*$v*(1-cos($theta)) + $w*sin($theta)]
    set M23 [expr $v*$w*(1-cos($theta)) - $u*sin($theta)]
    set M31 [expr $u*$w*(1-cos($theta)) - $v*sin($theta)]
    set M32 [expr $v*$w*(1-cos($theta)) + $u*sin($theta)]

    set M14 [expr (($a*($v**2 + $w**2) - $u*($b*$v + $c*$w)) * (1-cos($theta))) + (($b*$w - $c*$v) * sin($theta))]
    set M24 [expr (($b*($u**2 + $w**2) - $v*($a*$u + $c*$w)) * (1-cos($theta))) + (($c*$u - $a*$w) * sin($theta))]
    set M34 [expr (($c*($u**2 + $v**2) - $w*($a*$u + $b*$v)) * (1-cos($theta))) + (($a*$v - $b*$u) * sin($theta))]

    set row1 [list $M11 $M12 $M13 $M14]
    set row2 [list $M21 $M22 $M23 $M24]
    set row3 [list $M31 $M32 $M33 $M34]
    set row4 [list  0.0  0.0  0.0  1.0]

    set M [list $row1 $row2 $row3 $row4]

    return $M
}

# ============================================================================ #

proc rotate_phi {molID residue angle nf} {
    # rotate one phi dihedral angle of frame nf
    set pointSel  [atomselect $molID "residue $residue and name N" frame $nf]
    set point     [lindex [$pointSel get {x y z}] 0]
    set vecSel    [atomselect $molID "residue $residue and name CA" frame $nf]
    set direction [vecsub [lindex [$vecSel get {x y z}] 0] $point]
#    mol modselect 0 $molID "residue $residue and name N CA"

#    puts "$point $direction"

    set M [rotation_matrix_arbitrary_line $point $direction $angle]

    #set rotSel    [atomselect $molID "(residue > $residue) or (residue $residue and not name N H CA(sidechain or name HA C O))" frame $nf]
    set rotSel    [atomselect $molID "(residue > $residue) or (residue $residue and not name N H CA)" frame $nf]
#    mol modselect 4 $molID [$rotSel text]
    $rotSel move $M
}

# ============================================================================ #

proc rotate_psi {molID residue angle nf} {
    # rotate one psi dihedral angle of frame nf
    set pointSel  [atomselect $molID "residue $residue and name CA" frame $nf]
    set point     [lindex [$pointSel get {x y z}] 0]
    set vecSel    [atomselect $molID "residue $residue and name C" frame $nf]
    set direction [vecsub [lindex [$vecSel get {x y z}] 0] $point]

    set M [rotation_matrix_arbitrary_line $point $direction $angle]

    #set rotSel    [atomselect $molID "(residue > $residue) or (residue $residue and name O)" frame $nf]
    set rotSel    [atomselect $molID "(residue > $residue) or (residue $residue and not (sidechain or name N H H2 H3 CA HA C))" frame $nf]
    $rotSel move $M
}

# ============================================================================ #

proc rotate_dihedrals {molID rotation_vec f nf} {
    # rotate all phi and psi angles of frame nf by fraction f of rotation_vec

    lassign $rotation_vec phi_vec psi_vec
    set phi_vec [vecscale $f $phi_vec]
    set psi_vec [vecscale $f $psi_vec]

    set CAsel [atomselect $molID "name CA" frame $nf]
    set nres [llength [$CAsel get residue]]
#    set nres 2
    for {set i 0} {$i < $nres} {incr i} {
#    for {set i 99} {$i < $nres} {incr i} {}
#    set residue 2
#    for {set i $residue} {$i < [expr $residue + 1]} {incr i} {}
        if {$i > 0} {
            rotate_phi $molID $i [lindex $phi_vec [expr $i - 1]] $nf
        }
        if {$i < [expr $nres - 1]} {
            rotate_psi $molID $i [lindex $psi_vec $i] $nf
        }
    }
}

# ============================================================================ #

proc interpolate_sidechains {molID f nf} {
    # interpolate sidechain coordinates
    set CAsel   [atomselect $molID "name CA" frame $nf]
    set nres    [llength [$CAsel get residue]] 
    set nframes [molinfo $molID get numframes]

    for {set residue 0} {$residue < $nres} {incr residue} {
        # select sidechain of first, last and current frames
        set selFirstCA [atomselect $molID "residue $residue and (sidechain or name CA HA)" frame 0]
        set selLastCA  [atomselect $molID "residue $residue and (sidechain or name CA HA)" frame [expr $nframes - 1]]
        set selNowCA   [atomselect $molID "residue $residue and (sidechain or name CA HA)" frame $nf]
        set selFirst   [atomselect $molID "residue $residue and (sidechain or name CA HA)" frame 0]
        set selLast    [atomselect $molID "residue $residue and (sidechain or name CA HA)" frame [expr $nframes - 1]]
        set selNow     [atomselect $molID "residue $residue and (sidechain or name CA HA)" frame $nf]

        # backup first frame coordinates
        set firstFrameBackup [$selFirst get {x y z}]

        # superimpose first frame on last
        set M [measure fit $selFirstCA $selLastCA]
        $selFirst move $M
        
        # compute scaled difference vector
        set xFirst [$selFirst get x]
        set yFirst [$selFirst get y]
        set zFirst [$selFirst get z]
        set xLast  [$selLast  get x]
        set yLast  [$selLast  get y]
        set zLast  [$selLast  get z] 
#        set xdiff  [vecscale $f [vecsub $xLast $xFirst]]
#        set ydiff  [vecscale $f [vecsub $yLast $yFirst]]
#        set zdiff  [vecscale $f [vecsub $zLast $zFirst]]

        # interplate first frame coordinates
        $selFirst set x [vecadd [vecscale [expr {1.0 - $f}] $xFirst] [vecscale $f $xLast]]
        $selFirst set y [vecadd [vecscale [expr {1.0 - $f}] $yFirst] [vecscale $f $yLast]]
        $selFirst set z [vecadd [vecscale [expr {1.0 - $f}] $zFirst] [vecscale $f $zLast]] 

        # fit first frame to current
        set M [measure fit $selFirstCA $selNowCA]
        $selFirst move $M

        # transfer coordinates
        $selNow set x [$selFirst get x]
        $selNow set y [$selFirst get y]
        $selNow set z [$selFirst get z]

        # update current frame coordinates
#        set xNow [$selNow get x]
#        set yNow [$selNow get y]
#        set zNow [$selNow get z] 
#        $selNow set x [vecadd $xFirst $xdiff]
#        $selNow set y [vecadd $yFirst $ydiff]
#        $selNow set z [vecadd $zFirst $zdiff]
        
        # restore first frame coordinates
        $selFirst set {x y z} $firstFrameBackup
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


    # Make N-2 new frames (copied from the last frame)
    for {set i 2} {$i < $N} {incr i} {
	    animate dup frame 1 $molid
    }
    # there are now N frames
    # </plagiarism>
    
    set sumf 0
   
    # <plagiarism>
    # Do the linear interpolation in steps of 1/N so
    # f(0) = 0.0 and f(N-1) = 1.0
    for {set t 1} {$t < [expr $N - 1]} {incr t} {

        puts -nonewline [format "\rProgress: %5.1f%%" [expr 100.0*$t/($N-2)]]
        flush stdout

	    # Here's the call to the user-defined morph function
	    set f [$morph_type $t [expr $N - 2]]; # f grows from 0 to 1
        # </plagiarism>
	    set fprev [$morph_type [expr $t - 1] [expr $N - 2]]; # f grows from 0 to 1
        set df   [expr $f - $fprev]
        set sumf [expr $sumf + $df]
        set f [expr 1.0 / ($N - 1 - $t)]

        # copy coordinates from previous frame
        set currentFrame  [atomselect $molid "all" frame $t]
        set previousFrame [atomselect $molid "all" frame [expr $t - 1]]
        $currentFrame set {x y z} [$previousFrame get {x y z}]

        # interpolate dihedrals
        set rotation_vec [get_rotation_vector $molid $t [expr $N - 1]]
        lassign $rotation_vec phi_vec psi_vec
        set phiminmax [vecminmax $phi_vec]
        set psiminmax [vecminmax $psi_vec]
        puts "[format " distance: %16.12f, df: %6.4f, sumf: %6.4f, f: %5.3f, minmax: %5.1f, %5.1f" [veclength $phi_vec] $df $sumf $f [lindex $phiminmax 0] [lindex $phiminmax 1]]"

        # do the rotations
        rotate_dihedrals $molid $rotation_vec $f $t

        # interpolate sidechain coordinates
        interpolate_sidechains $molid $f $t

        # fit current frame to previous frame
        set previousFrame [atomselect $molid "all" frame [expr $t - 1]]
        set M [measure fit $currentFrame $previousFrame]
        $currentFrame move $M

#        # compute distance from target
#        set current_rotation_vec [get_rotation_vector $molid $t [expr $N - 1]]
#        lassign $current_rotation_vec current_phi_vec current_psi_vec
#        puts $current_phi_vec
        #puts [vecmean $current_phi_vec]
        #puts "[vecminmax $current_phi_vec] [vecminmax $current_psi_vec]"
    }
    puts ""

    # final rotations
    set current_rotation_vec [get_rotation_vector $molid [expr $N - 2] [expr $N - 1]]
    rotate_dihedrals $molid $current_rotation_vec 1.0 [expr $N - 2]

    # fit current frame to previous frame
    set currentFrame  [atomselect $molid "all" frame [expr $N - 2]]
    set previousFrame [atomselect $molid "all" frame [expr $N - 3]]
    set M [measure fit $currentFrame $previousFrame]
    $currentFrame move $M 

    set current_rotation_vec [get_rotation_vector $molid [expr $N - 2] [expr $N - 1]]
    lassign $current_rotation_vec current_phi_vec current_psi_vec
    puts $current_phi_vec 
    puts ""
    puts $current_psi_vec 

    # fit final frame to final but one frame
    set finalFrame    [atomselect $molid "all" frame last]
    set previousFrame [atomselect $molid "all" frame [expr $N - 2]]
    set M [measure fit $finalFrame $previousFrame]
    $finalFrame move $M 


    animate goto 0



#    puts $current_phi_vec
}


proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.15
    graphics $mol cone $middle $end radius 0.25
}


proc vecminmax {v} {
    set max [lindex $v 0]
    set min [lindex $v 0]
    foreach n $v {
        if {$n > $max} {
            set max $n
        }
        if {$n < $min} {
            set min $n
        }
    }
    return [list $min $max]
}
