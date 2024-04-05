#!/usr/bin/env perl

# This script goes through the DST file output by DynaNet after performing a
# near station search. It removes stations that jurisdictions have listed in
# in their .near files

# Read in the APREF stations
for (glob'/home/fedora/apref*.disconts') {$file = $_};
open APREF, $file;
while (<APREF>) {
    unless (/^\#/) {
        @_ = split ' ';
        $aprefStation{$_[0]}++;
    }
}

# Read in the true near stations
for $file (glob'*.near') {
    open IN, $file;
    while (<IN>) {
        $stat1 = substr $_, 0, 20;
        $stat2 = substr $_, 20, 20;
        $stat1 =~ s/^\s+//;
        $stat1 =~ s/\s+$//;
        $stat2 =~ s/^\s+//;
        $stat2 =~ s/\s+$//;
        unless ($aprefStation{$stat1} && $aprefStation{$stat2}) {    
            $ignore{"$stat1 $stat2"}++;
            $correct[$numIgnore][0] = $stat1;
            $correct[$numIgnore][1] = $stat2;
            $numIgnore++;
        }
    }
}
print "There are $numIgnore pairs of near stations to ignore\n";

# Find the .dst file
for $dstFile (glob'*.dst') {
    push @dstFiles, $dstFile;
}
$numFiles = $#dstFiles + 1;

# If there is multiple .dst files allow the user to determine which one to use
# Die if there is no .dst file
if ($#dstFiles == 0) {
    $file = $dstFiles[0];
} elsif ($#dstFiles < 0) {
    die "There is no dst file to check.\n";
} else {
    print "There are multiple dst files:\n";
    for (@dstFiles) {
        $i++;
        print "\t$i.\t$_\n";
    }
    print "Type the number of the file you want to check: ";
    chomp($fileNum = <STDIN>);
    if ($fileNum < 1 || $fileNum > $numFiles) {
        print 'Invalid response. Please select a number between 1 and ';
        print "$numFiles\n";
        die "\n";
    }
    $file = $dstFiles[$fileNum - 1];
}

# loop over the .dst file and write out pairs that aren't included in any .near
# file
$bakFile = $file . '.bak';
`mv $file $bakFile`;
open OUT, ">$file";
open IN, $bakFile;
while (<IN>) {
    last if (/^First station/);
}
$_ = <IN>;
while (<IN>) {
    unless (/^\n/) {
        $numNear++;
        $line = $_;
        $stat1 = substr $_, 0, 20;
        $stat2 = substr $_, 20, 20;
        $stat1 =~ s/^\s+//;
        $stat1 =~ s/\s+$//;
        $stat2 =~ s/^\s+//;
        $stat2 =~ s/\s+$//;
        $printIt = 1;
        if ($aprefStation{$stat1} && $aprefStation{$stat2}) {
            $numAPREF++;
        } else {
            for (0..$numIgnore-1) {
                if ($stat1 eq $correct[$_][0] && $stat2 eq $correct[$_][1]) {
                    $printIt = 0;
                    $ignore{"$stat1 $stat2"}--;
                    $numMatches++;
                    last;
                }
                if ($stat1 eq $correct[$_][1] && $stat2 eq $correct[$_][0]) {
                    $printIt = 0;
                    $ignore{"$stat2 $stat1"}--;
                    $numMatches++;
                    last;
                } 
            }
            if ($printIt) {
                print OUT $line;
                $numLeft++;
            }
        }            
    }
}
close IN;
close OUT;
print "There are $numNear pairs of near stations\n";
print "There are $numMatches of these near station pairs in the ignore files\n";
print "There should be " . ($numNear - $numMatches - $numAPREF) .
    " near matches ";
print "left and...";
if (($numNear - $numMatches - $numAPREF) == $numLeft) {
    print "there is!\n";
} else {
    print "there's not. :-/ There's $numLeft left\n";
}
open OUT, '>notUsedIgnore.dat'; 
for $key (keys %ignore) {
    if ($ignore{$key} == 1) {
        $ignoreLeft++;
        print OUT "$key\n";
    }
}
print "There should be " . ($numIgnore - $numMatches) . " near matches still ";
print "to ignore and...";
if (($numIgnore - $numMatches) == $ignoreLeft) {
    print "there is!\n";
} else {
    print "there's not. :-/ There's $ignoreLeft left\n";
}
