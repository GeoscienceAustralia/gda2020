#!/usr/bin/env perl

# This script runs a minimally-constrained adjustment for each baseline cluster

# If necessary, copy over the files needed to do the adjustment
`cp ~/DynaML.xsd .` unless (-e 'DynaML.xsd');
unlink glob'apref*.snx';
`cp ~/apref/apref*.snx .`;
unlink glob'disconts*.snx';
`cp ~/apref/disconts*.snx .`;
for (glob'disconts*.snx') {$discontsFile = $_}

# Create a list of APREF stations
for (glob'apref*.snx') {$file = $_}
@apref_stns = ();
open APREF, $file;
while (<APREF>) {
    if (/-SITE\/ID/) {last}
    if ($go) {push @apref_stns, substr($_, 1, 4)}
    if (/\+SITE\/ID/) {$go++}
}

# Loop over the XML baseline files
for $stnFile (glob'*_all_stn.xml') {
    ($msrFile = $stnFile) =~ s/stn/msr/;
    ($network = $stnFile) =~ s/\_stn\.xml//;
    %done = ();
    $excl = '';
    $incl = '--include-stns-assoc-msrs "';
    $constraints = '--constraints "';
    open STN, $stnFile;
    $mac1 = 0;
    while (<STN>) {
    	if (/\<Name\>/) {
            ($stn = $_) =~ s/\s+\<Name\>//;
            $stn =~ s/\<\/Name\>\s+//;
            if ($stn =~ /^MAC1\_\d{7}$/) {
                $mac1 = 1;
                $excl = '--exclude-stns-assoc-msrs "' . $stn . '"';
            }
            unless ($done{$stn}++ or $stn =~ /^MAC1/) {
            	$incl .= "$stn,";
                if (grep(/^$stn$/,@apref_stns)) {
                    $constraints .= "$stn,CCC,";
                } elsif ($stn =~ /\w{4}\_\d{7}/) {
                    $root = substr($stn, 0, 4);
                    if (grep(/^$root$/,@apref_stns)) {
                        $constraints .= "$stn,CCC,";
                    }
                }
            }
        }
    }
    open(FH, '>', $network);   
    $incl =~ s/,$/"/;
    $constraints =~ s/,$/"/;
    $cmd = "dnaimport -n $network -r GDA2020 apref*.snx $stnFile ";
    $cmd = $cmd . "$msrFile $incl --split-gnss-cluster-msrs ";
    if ($mac1 == 1) {
        $cmd = $cmd . "$excl ";
    }
    $cmd = $cmd . "--exclude-msr-types 'Y' --discontinuity-file $discontsFile";
    print FH "$cmd \n";
    `$cmd`;    
    $cmd = '';
    $cmd = "dnareftran -n $network -r GDA2020";
    print FH "$cmd \n";
    `$cmd`;
    $cmd = '';
    $cmd = "dnaadjust -n $network --simult $constraints ";
    $cmd = $cmd . "--output-adj-msr --sort-adj-msr-field 7";
    print FH "$cmd \n";
    `$cmd`;
}
