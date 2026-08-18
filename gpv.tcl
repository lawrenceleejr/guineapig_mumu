#!/usr/local/bin/wish

set fname [lindex $argv 0]
set lh {}
foreach n [exec gpv $fname] {
    lappend lh $n
}
set cmd [concat exec gpv $fname $lh]
puts [eval $cmd]
exit
