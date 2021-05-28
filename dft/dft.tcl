proc dft {x} {
    # TCL doesn't support complex numbers so will keep them separate and use Euler's formula
    set N [llength $x]
    set pi [expr "acos(-1)"]
    set Hk {}
    for {set k 0} {$k < $N} {incr k} {
        set sum_real 0.0
        set sum_imag 0.0
        for {set n 0} {$n < $N} {incr n} {
            set xn [lindex $x $n]
            set real [expr {$xn*cos(-2*$pi*$n*$k/double($N))}]
            set imag [expr {$xn*sin(-2*$pi*$n*$k/double($N))}]
            set sum_real [expr {$sum_real + $real}]
            set sum_imag [expr {$sum_imag + $imag}]
        }
        lappend Hk $sum_real $sum_imag
    }
    return $Hk
}

proc listFromFile {filename} {
    set f [open $filename r]
    set data [split [string trim [read $f]]]
    close $f
    return $data
}

set x [listFromFile data.csv]
set N [llength $x]
puts "Read in $N lines for vsin"

puts "Calculating DFT"
set t0 [clock milliseconds]
set Hk [dft $x]
set t1 [clock milliseconds]
set dt [expr {($t1-$t0)/1000.0}]
puts "$N-point DFT took $dt seconds."

proc print_Hk {Hk} {
    foreach {real imag} $Hk {
        puts "$real $imag"
    }
}