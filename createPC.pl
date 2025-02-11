#!/usr/bin/env perl

# The script creates a single DynaML GNSS point cluster (Type Y) from the
# times series combination SINEX file. It also sets Vscale to 40, which was 
# empirically determined by a comparison of the APREF weekly solution standard
# deviations to those given in the Reg13 certificates

use DateTime::Precise;

# Create dictionaries to hold the discontinuity epoch data. These data are read
# in from the SOLUTION/DISCONTINUITY block in soln.snx. This file holds data
# on all IGS stations, not just Australian stations
open DISCONT, 'soln.snx';
while (<DISCONT>) {
    $discontLine = $_;
    last if ($discontLine =~ /\-SOLUTION\/DISCONTINUITY/);
    if ($go) {
        @_ = split ' ', $_;
        if ($_[6] eq 'P') {
            $segment = $_[0] . $_[2];
            $allYearDOY{$segment} = substr($discontLine, 16, 2) . 
                substr($discontLine, 19, 3);
            if ($allYearDOY{$segment} < 94000) {
                $allYearDOY{$segment} = '20' . $allYearDOY{$segment};
            } else {
                $allYearDOY{$segment} = '19' . $allYearDOY{$segment};
            }
        }                
    }
    $go = 1 if ($discontLine =~ /\+SOLUTION\/DISCONTINUITY/);
}
$go = 0;

# Open and read in SNXEPO.SNX
open SNX, 'SNXEPO.SNX';
for (<SNX>) {push @snxLines, $_}

# Open the output files
$aprefEpoch = substr $snxLines[0], 32, 6;
$aprefEpoch =~ s/://;
$stn = 'apref' . $aprefEpoch . '_stn.xml';
$msr = 'apref' . $aprefEpoch . '_msr.xml';
open STN, "> $stn";
open MSR, "> $msr";

# Get the solution epoch data. These data will be used to remove non-Australian
# stations from the discontinuity epoch data 
for $snxLine (@snxLines) {
    last if ($snxLine =~ /\-SOLUTION\/EPOCHS/);
    if ($go && !($snxLine =~ /^\*Code/)) {
        push @epochs, $snxLine;
    }
    $go = 1 if ($snxLine =~ /\+SOLUTION\/EPOCHS/);
}
$go = 0;

# Get the station estimates
for $snxLine (@snxLines) {
    last if ($snxLine =~ /\-SOLUTION\/ESTIMATE/);
    if ($go && !($snxLine =~ /^\*INDEX/)) {
        push @estimates, $snxLine;
    }
    $go = 1 if ($snxLine =~ /\+SOLUTION\/ESTIMATE/);
}
$go = 0;

# Go through solution epoch data in order to filter out non-Australian stations
# from the discontinuity epoch data. Also, replace the epoch for the first
# point code, i.e. 2000000, with the starting epoch of the solution epoch data. 
# This is a quirk of catref 
for $epoch (@epochs) {
    @_ = split ' ', $epoch;
    $segment = $_[0] . $_[2];
    $stat{$segment} = $_[0];
    $pointCode{$segment} = $_[2];
    if ($allYearDOY{$segment}) {
        if ($allYearDOY{$segment} eq '2000000') {
            $yearDOY{$segment} = '1900001';
        } else {
            $yearDOY{$segment} = $allYearDOY{$segment}
        }
    } else {
        $yearDOY{$segment} = substr($epoch, 16, 2) . substr($epoch, 19, 3);
        if ($yearDOY{$segment} < 94000) {
            $yearDOY{$segment} = '20' . $yearDOY{$segment};
        } else {
            $yearDOY{$segment} = '19' . $yearDOY{$segment};
        }    
    }    
    $numSegs{$stat{$segment}}++;
}

# Create dictionaries to hold the estimate data
for $estimate (@estimates) {
    @_ = split ' ', $estimate;

# Get reference epoch
    unless ($gotEpoch++) {
        ($yr, $doy) = split(':', $_[5]);
        $doy--;
        $epoch = '20' . $yr . '.01.01'; 
        $dt = DateTime::Precise->new($epoch);
        $dt->inc_day($doy);
        $epoch = $dt->strftime('%d.%m.%Y');
    }
    $segment = $_[2] . $_[4];
    $stax{$segment} = $_[8] if ($_[1] eq 'STAX');
    $stddevx{$segment} = $_[9] if ($_[1] eq 'STAX');
    $stay{$segment} = $_[8] if ($_[1] eq 'STAY');
    $stddevy{$segment} = $_[9] if ($_[1] eq 'STAY');
    $staz{$segment} = $_[8] if ($_[1] eq 'STAZ');
    $stddevz{$segment} = $_[9] if ($_[1] eq 'STAZ');
}

# Output the information (sorted alphabetically) to aprefDisconts.dat
open TMP1, '>tmp1';
for $key (keys %stat) {
    # Below is because stations like KRNG has only one segement despite having
    # a discontinuity because not all the historicla data has been processed 
    # yet
    if ($numSegs{$stat{$key}} == 1 && substr($key, -1) eq '1') {
        $stat = $stat{$key};
        print TMP1 sprintf("%-15s %2s %21s %21s %21s %11s %11s %11s\n", $stat,
            $pointCode{$key}, $stax{$key}, $stay{$key}, $staz{$key},
            $stddevx{$key}, $stddevy{$key}, $stddevz{$key});
    } else {
        $stat = "$stat{$key}\_$yearDOY{$key}";
        print TMP1 sprintf("%-15s %2s %21s %21s %21s %11s %11s %11s\n", $stat,
            $pointCode{$key}, $stax{$key}, $stay{$key}, $staz{$key},
            $stddevx{$key}, $stddevy{$key}, $stddevz{$key});
    }
    $stn{$key} = $stat;
}
$out = 'apref'. $aprefEpoch . '.disconts';
open OUT, ">$out";
print OUT "#Stn           Pt X                     Y                    Z";
print OUT "                     SDX         SDY         SDZ\n";
close OUT;
`sort tmp1 >> $out`;
unlink 'tmp1';

# Write out the station data to file
print STN '<?xml version="1.0"?>' . "\n";
print STN '<DnaXmlFormat type="Station File" referenceframe="GDA2020" ';
print STN 'epoch="';
print STN $epoch;
print STN '" ';
print STN 'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ';
print STN 'xsi:noNamespaceSchemaLocation="DynaML.xsd">' . "\n";
@stn = @stax = @stay = @staz = ();
for $snxLine (@snxLines) {
    last if ($snxLine =~ /\-SOLUTION\/ESTIMATE/);
    if ($go && !($snxLine =~ /^\*INDEX/)) {
        @col = split ' ', $snxLine;
        $col[8] =~ s/^\-\./-0./;
        if ($col[1] eq 'STAX') {
            $segment = "$col[2]$col[4]";
            push @stn, $stn{$segment};
            push @stax, $col[8];
        }
        push @stay, $col[8] if ($col[1] eq 'STAY');
        push @staz, $col[8] if ($col[1] eq 'STAZ');
    }
    $go = 1 if ($snxLine =~ /\+SOLUTION\/ESTIMATE/);
}
$go = 0;
for $i (0..$#stn) {
    print STN "\t<DnaStation>\n";
    print STN "\t\t<Name>$stn[$i]</Name>\n";
    print STN "\t\t<Constraints>FFF</Constraints>\n";
    print STN "\t\t<Type>XYZ</Type>\n";
    print STN "\t\t<StationCoord>\n";
    print STN "\t\t\t<Name>$stn[$i]</Name>\n";
    print STN "\t\t\t<XAxis>$stax[$i]</XAxis>\n";
    print STN "\t\t\t<YAxis>$stay[$i]</YAxis>\n";
    print STN "\t\t\t<Height>$staz[$i]</Height>\n";
    print STN "\t\t</StationCoord>\n";
    print STN "\t\t<Description>$stn[$i]</Description>\n";
    print STN "\t</DnaStation>\n";
}
print STN "</DnaXmlFormat>\n";
close STN;

# Write out the measurement data to file
print MSR '<?xml version="1.0"?>' . "\n";
print MSR '<DnaXmlFormat type="Measurement File" referenceframe="GDA2020" ';
print MSR 'epoch="';
print MSR $epoch;
print MSR '" ';
print MSR 'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ';
print MSR 'xsi:noNamespaceSchemaLocation="DynaML.xsd">' . "\n";
for $snxLine (@snxLines) {
    last if ($snxLine =~ /\-SOLUTION\/MATRIX_ESTIMATE/);
    if ($go && !($snxLine =~ /^\*PARA1/)) {
        @cols = split ' ', $snxLine;
        $row = $cols[0] - 1;
        $col0 = $cols[1] - 1;
        for $i (0..$#cols-2) {
            $col = $col0 + $i;
            $cols[$i+2] =~ s/^\-\./-0./;
            $vcv[$row][$col] = $cols[$i+2];
        }
    }
    if ($snxLine =~ /\+SOLUTION\/MATRIX_ESTIMATE/) {
        $lower++ if ($snxLine =~ /L COVA/);
        $upper++ if ($snxLine =~ /U COVA/);
        $go = 1;
    }        
}
$go = 0;
$setNumber = $#stn + 1;
print MSR "\t<!--Type Y GPS Point Cluster (set of $setNumber)-->\n";
print MSR "\t<DnaMeasurement>\n";
print MSR "\t\t<Type>Y</Type>\n";
print MSR "\t\t<Ignore/>\n";
print MSR "\t\t<ReferenceFrame>GDA2020</ReferenceFrame>\n";
print MSR "\t\t<Epoch>$epoch</Epoch>\n";
print MSR "\t\t<Vscale>1.000</Vscale>\n";
print MSR "\t\t<Pscale>1.000</Pscale>\n";
print MSR "\t\t<Lscale>1.000</Lscale>\n";
print MSR "\t\t<Hscale>1.000</Hscale>\n";
print MSR "\t\t<Coords>XYZ</Coords>\n";
print MSR "\t\t<Total>$setNumber</Total>\n";
if ($lower) {
    for $i (0..$#stn) {
        print MSR "\t\t<First>$stn[$i]</First>\n";
        print MSR "\t\t<Clusterpoint>\n";
        print MSR "\t\t\t<X>$stax[$i]</X>\n";
        print MSR "\t\t\t<Y>$stay[$i]</Y>\n";
        print MSR "\t\t\t<Z>$staz[$i]</Z>\n";
        print MSR "\t\t\t<SigmaXX>$vcv[6*$i][6*$i]</SigmaXX>\n";
        print MSR "\t\t\t<SigmaXY>$vcv[6*$i+1][6*$i]</SigmaXY>\n";
        print MSR "\t\t\t<SigmaXZ>$vcv[6*$i+2][6*$i]</SigmaXZ>\n";
        print MSR "\t\t\t<SigmaYY>$vcv[6*$i+1][6*$i+1]</SigmaYY>\n";
        print MSR "\t\t\t<SigmaYZ>$vcv[6*$i+2][6*$i+1]</SigmaYZ>\n";
        print MSR "\t\t\t<SigmaZZ>$vcv[6*$i+2][6*$i+2]</SigmaZZ>\n";
        for $j ($i+1..$#stn) {
            print MSR "\t\t\t<PointCovariance>\n";
            print MSR "\t\t\t\t<m11>$vcv[6*$j][6*$i]</m11>\n";
            print MSR "\t\t\t\t<m12>$vcv[6*$j+1][6*$i]</m12>\n";
            print MSR "\t\t\t\t<m13>$vcv[6*$j+2][6*$i]</m13>\n";
            print MSR "\t\t\t\t<m21>$vcv[6*$j][6*$i+1]</m21>\n";
            print MSR "\t\t\t\t<m22>$vcv[6*$j+1][6*$i+1]</m22>\n";
            print MSR "\t\t\t\t<m23>$vcv[6*$j+2][6*$i+1]</m23>\n";
            print MSR "\t\t\t\t<m31>$vcv[6*$j][6*$i+2]</m31>\n";
            print MSR "\t\t\t\t<m32>$vcv[6*$j+1][6*$i+2]</m32>\n";
            print MSR "\t\t\t\t<m33>$vcv[6*$j+2][6*$i+2]</m33>\n";
            print MSR "\t\t\t</PointCovariance>\n";
        }
        print MSR "\t\t</Clusterpoint>\n";
    }
} elsif ($upper) {
    for $i (0..$#stn) {
        print MSR "\t\t<First>$stn[$i]</First>\n";
        print MSR "\t\t<Clusterpoint>\n";
        print MSR "\t\t\t<X>$stax[$i]</X>\n";
        print MSR "\t\t\t<Y>$stay[$i]</Y>\n";
        print MSR "\t\t\t<Z>$staz[$i]</Z>\n";
        print MSR "\t\t\t<SigmaXX>$vcv[6*$i][6*$i]</SigmaXX>\n";
        print MSR "\t\t\t<SigmaXY>$vcv[6*$i][6*$i+1]</SigmaXY>\n";
        print MSR "\t\t\t<SigmaXZ>$vcv[6*$i][6*$i+2]</SigmaXZ>\n";
        print MSR "\t\t\t<SigmaYY>$vcv[6*$i+1][6*$i+1]</SigmaYY>\n";
        print MSR "\t\t\t<SigmaYZ>$vcv[6*$i+1][6*$i+2]</SigmaYZ>\n";
        print MSR "\t\t\t<SigmaZZ>$vcv[6*$i+2][6*$i+2]</SigmaZZ>\n";
        for $j ($i+1..$#stn) {
            print MSR "\t\t\t<PointCovariance>\n";
            print MSR "\t\t\t\t<m11>$vcv[6*$i][6*$j]</m11>\n";
            print MSR "\t\t\t\t<m12>$vcv[6*$i][6*$j+1]</m12>\n";
            print MSR "\t\t\t\t<m13>$vcv[6*$i][6*$j+2]</m13>\n";
            print MSR "\t\t\t\t<m21>$vcv[6*$i+1][6*$j]</m21>\n";
            print MSR "\t\t\t\t<m22>$vcv[6*$i+1][6*$j+1]</m22>\n";
            print MSR "\t\t\t\t<m23>$vcv[6*$i+1][6*$j+2]</m23>\n";
            print MSR "\t\t\t\t<m31>$vcv[6*$i+2][6*$j]</m31>\n";
            print MSR "\t\t\t\t<m32>$vcv[6*$i+2][6*$j+1]</m32>\n";
            print MSR "\t\t\t\t<m33>$vcv[6*$i+2][6*$j+2]</m33>\n";
            print MSR "\t\t\t</PointCovariance>\n";
        }
        print MSR "\t\t</Clusterpoint>\n";
    }
} else {
    die "VCV is not a lower nor an upper triangular matrix\n";
}    
print MSR "\t</DnaMeasurement>\n";
print MSR "</DnaXmlFormat>\n";
close MSR;
