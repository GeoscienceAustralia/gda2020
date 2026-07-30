#!/usr/bin/env perl

# This script will examine the stations listed in a dst file to pull out
# genuine duplicates

for $dstFile (glob'*.dst') {
    push @dstFiles, $dstFile;
}
$numFiles = $#dstFiles + 1;

if ($#dstFiles == 0) {
    $dFile = $dstFiles[0];
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
    $dFile = $dstFiles[$fileNum - 1];
}

for (glob'~/apref*.disconts') {$file = $_};
open APREF, $file;
while (<APREF>) {
    unless (/^\#/) {
	@_ = split ' ';
        $aprefStation{$_[0]}++;
    }
}

for $file (glob'~/renaming/*.renaming') {
    open IN, $file;
    while (<IN>) {
        unless (/^!\#/) {
            $station = substr $_, 0, 20;
            $station =~ s/^\s*//;
            $station =~ s/\s*$//;
            $renameStation{$station}++;
        }
    }
}

for $file (glob'~/renaming/*.ignore') {
    open IN, $file;
    while (<IN>) {
        chomp;
        s/^\s*//;
        s/\s*$//;
        $ignoreStation{$_}++;
    }
}

for $stn (glob'stn/inApriori/*.xml') {
    open IN, $stn;
    @_ = split '_', $stn;
    ($juris = $_[0]) =~ s/^stn\/inApriori\///;
    while (<IN>) {
        if (/<DnaStation>/) {
            chomp($station = <IN>);
            $station =~ s/^\s*<Name>//;
            $station =~ s/<\/Name>\s*$//;
            push @{$statJuris{$station}}, $juris;
        }
    }
}

for $stn (glob'stn/*.xml') {
    open IN, $stn;
    @_ = split '_', $stn;
    ($juris = $_[0]) =~ s/^stn\///;
    while (<IN>) {
        if (/<DnaStation>/) {
            chomp($station = <IN>);
            $station =~ s/^\s*<Name>//;
            $station =~ s/<\/Name>\s*$//;
            push @{$statJuris{$station}}, $juris;
        }
    }
}
open IN, $dFile;
while (<IN>) {last if (/^Default reference frame/)}
for (0..2) {<IN>}
while (<IN>) {
    unless (/^\n/) {
        chomp;
        $_ = substr $_, 4;
        unless ($aprefStation{$_} or $renameStation{$_} or $ignoreStation{$_}) {
            push @dups, $_;
        }
    }
}


if (@dups) {
    ($out = $dFile) =~ s/dst/dup/;
    open OUT, ">$out";
    print "Duplicates were found and written to $out\n";
    for (@dups) {
        print OUT "$_";
        foreach (@{$statJuris{$_}}) {
            print OUT "\t$_";
        }
        print OUT "\n";
    }
} else {
    print "Too sweet! No duplicates\n";
}

