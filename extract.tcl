#!/usr/local/bin/tclsh

proc getlumi {name} {
    set file [open $name r]
    gets $file line
    set limit0 [expr 0.0*3000.0]
    set sum0 0.0
    set limit1 [expr 0.95*3000.0]
    set sum1 0.0
    set limit2 [expr 0.99*3000.0]
    set sum2 0.0
    while {[eof $file ]==0} {
	gets $file lineold
	gets $file line
	if {[eof $file]==0} {
	    if {[lindex $lineold 0]>=$limit0} {
		set sum0 [expr $sum0+[lindex $line 1]]
		if {[lindex $lineold 0]>=$limit1} {
		    set sum1 [expr $sum1+[lindex $line 1]]
		    if {[lindex $lineold 0]>=$limit2} {
			set sum2 [expr $sum2+[lindex $line 1]]
		    }
		}
	    }
	}
    }
    close $file
    return "$sum0 $sum1 $sum2"
}

set name clic.betax
#set name clic.offset_x

for {set i 6} {$i<20} {incr i} {
    set cmd [concat exec gpv $name.$i lumi_ee]
    eval $cmd
    puts "[expr $i*0.1] [getlumi lumi_ee.dat]"
}
