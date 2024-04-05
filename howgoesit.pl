#!/usr/bin/env perl

# This script reports on the progress of minAdjust.pl

$str1 = `ls *all_msr.xml | wc -l`;
$str2 = `ls *.adj | wc -l`;
chomp($str1);
chomp($str2);

$str3 = $str1 - $str2;

if ($str3 == 0) {
    print "All done!\n";
    `dynaclean.pl`;
    `grep Rigorous *.adj | sort -g -k 4 > sigma0.dat`;
    `grep Degrees *.adj | sort -g -k 4 > dof.dat`;
} else {
    print "$str2 of $str1 clusters have been processed; there are $str3 to go!\n";
}
