#!/usr/bin/env perl

# This script goes through a DynaML station file and a measurement file and
# translates stations that use the conventional 4-character ID to the name
# commonly used by the state or jurisdiction.

# Read in the translation table
$d = `pwd`;
@d = split('/', $d);
$jur = $d[-2];
$jur = 'qld' if ($jur eq 'surat');
$camp = '/home/fedora/ngca/' . $jur . '/baselines/';
chdir '/home/fedora/transTables';
$transT = $jur . 'TransTable*.csv';
`cp $transT $camp`;
for (glob"$transT") {
    open TT, $_;
    while (<TT>) {
        unless (/^\#/) {
            s/\s+$//;
            @_ = split ',', $_;
            $nameTrans{$_[0]} = $_[1];
        }
    }
}

# Read in the APREF stations and their discontinuities
for (glob'/home/fedora/apref/disconts????????.snx') {
    open APREF, $_;
    while (<APREF>) {
        unless (/\*-----|%ENDSNX|-SOLUTION|\+SOLUTION/) {
            @_ = split ' ', $_;
            $apref{$_[0]}++;
        }
    }
}

# Loop over the files
chdir $camp;
for $stn (glob'*stn.xml') {
    ($msr = $stn) =~ s/stn/msr/;
    ($clus = $stn) =~ s/\_stn\.xml//;

# Make copies of the DynaML files
    `mv $stn tmpSTN`;
    `mv $msr tmpMSR`;

# Open the original DynaML station file and loop through it
    open OUT, ">$stn";
    open STN, 'tmpSTN';
    while (<STN>) {
        $line = $_;
        ($station = $line) =~ s/\s+//g;
        if ($line =~ /<Name>/) {
            $station =~ s/<Name>//;
            $station =~ s/<\/Name>//;
            unless ($apref{$station}) {
                $line =~ s/$station/$nameTrans{$station}/ 
                    if ($nameTrans{$station});
            }
        }
        if ($line =~ /<Description>/) {
            $station =~ s/<Description>//;
            $station =~ s/<\/Description>//;
            unless ($apref{$station}) {
                $line =~ s/$station/$nameTrans{$station}/ 
                    if ($nameTrans{$station});
            }
        }
        print OUT $line;
    }

# Open the original DynaML measurement file and loop through it
    open OUT, ">$msr";
    open MSR, 'tmpMSR';
    while (<MSR>) {
        $line = $_;
        ($station = $line) =~ s/\s+//g;
        if ($line =~ /<First>/) {
            $station =~ s/<First>//;
            $station =~ s/<\/First>//;
            unless ($apref{$station}) {
                $line =~ s/$station/$nameTrans{$station}/ 
                    if ($nameTrans{$station});
            }
        }
        if ($line =~ /<Second>/) {
            $station =~ s/<Second>//;
            $station =~ s/<\/Second>//;
            unless ($apref{$station}) {
                $line =~ s/$station/$nameTrans{$station}/ 
                    if ($nameTrans{$station});
            }
        }
        print OUT $line;
    }
}
unlink 'tmpSTN', 'tmpMSR';
