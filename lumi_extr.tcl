#!/usr/local/bin/tclsh

catch {exec rm lumi.ee}
catch {exec rm lumi.ep}
catch {exec rm lumi.pe}
catch {exec rm lumi.pp}
catch {exec rm lumi.eg}
catch {exec rm lumi.ge}
catch {exec rm lumi.pg}
catch {exec rm lumi.gp}
catch {exec rm lumi.gg}

for {set j 0} {$j<2} {incr j} {

    exec gunzip lumi.ee.$j.gz
    exec cp lumi.ee.$j lumi.ee.out
    exec gzip -9 lumi.ee.$j
    
    exec gunzip lumi.eg.$j.gz
    exec cp lumi.eg.$j lumi.eg.out
    exec gzip -9 lumi.eg.$j
    
    exec gunzip lumi.ge.$j.gz
    exec cp lumi.ge.$j lumi.ge.out
    exec gzip -9 lumi.ge.$j
    
    exec gunzip lumi.gg.$j.gz
    exec cp lumi.gg.$j lumi.gg.out
    exec gzip -9 lumi.gg.$j
    
    set file [open "lumi.ee.out" "r"]
    set f_ee [open "lumi.ee" "a"]
    set f_ep [open "lumi.ep" "a"]
    set f_pe [open "lumi.pe" "a"]
    set f_pp [open "lumi.pp" "a"]
    
    gets $file line
    while {![eof $file]} {
	set e1 [lindex $line 0]
	set e2 [lindex $line 1]
	set z [lindex $line 4]
	if {$e1<0.0} {
	    set e1 [string range $e1 1 end]
	    if {$e2<0.0} {
		set e2 [string range $e2 1 end]
		puts $f_pe "$e1 $e2 $z"
	    } {
		puts $f_pp "$e1 $e2 $z"
	    }
	} {
	    if {$e2<0.0} {
		set e2 [string range $e2 1 end]
		puts $f_ee "$e1 $e2 $z"
	    } {
		puts $f_ep "$e1 $e2 $z"
	    }
	}
	gets $file line
    }
    
    close $file
    close $f_ee
    close $f_ep
    close $f_pe
    close $f_pp
    
    set file [open "lumi.eg.out" "r"]
    set f_eg [open "lumi.eg" "a"]
    set f_pg [open "lumi.pg" "a"]
    
    gets $file line
    while {![eof $file]} {
	set e1 [lindex $line 0]
	set e2 [lindex $line 1]
	set z [lindex $line 4]
	if {$e1<0.0} {
	    set e1 [string range $e1 1 end]
	    puts $f_pg "$e1 $e2 $z"
	} {
	    puts $f_eg "$e1 $e2 $z"
	}
	gets $file line
    }
    
    close $file
    close $f_eg
    close $f_pg
    
    set file [open "lumi.ge.out" "r"]
    set f_ge [open "lumi.ge" "a"]
    set f_gp [open "lumi.gp" "a"]

    gets $file line
    while {![eof $file]} {
	set e1 [lindex $line 0]
	set e2 [lindex $line 1]
	set z [lindex $line 4]
	if {$e2<0.0} {
	    set e2 [string range $e2 1 end]
	    puts $f_ge "$e1 $e2 $z"
	} {
	    puts $f_gp "$e1 $e2 $z"
	}
	gets $file line
    }
    
    close $file
    close $f_ge
    close $f_gp
    
    set file [open "lumi.gg.out" "r"]
    set f_gg [open "lumi.gg" "a"]
    
    gets $file line
    while {![eof $file]} {
	set e1 [lindex $line 0]
	set e2 [lindex $line 1]
	set z [lindex $line 4]
	puts $f_gg "$e1 $e2 $z"
	gets $file line
    }
    
    close $file
    close $f_gg
}

exec rm lumi.ee.out
exec rm lumi.eg.out
exec rm lumi.ge.out
exec rm lumi.gg.out
