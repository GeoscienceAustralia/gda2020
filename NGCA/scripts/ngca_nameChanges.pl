#!/usr/bin/env perl

# This script changes the files that were renamed during processing back to
# their original name. Station names are also changed back in the SINEX files
# and RINEX file names in the antenna information files
# Run script from the archive directory

# Change renamed stations and RINEX files back to their original names
if (-e 'nameChanges.dat') {
    open IN, 'nameChanges.dat';
    while (<IN>) {
        @_ = split ' ', $_;
        push @newFile, $_[0];
        push @oldFile, $_[1];
        $rename{$_[0]} = $_[1];
    }    
    close IN;
    while (@newFile) {
        $newFile = pop @newFile;
        $newStat = substr $newFile, 0, 4;
        $oldFile = pop @oldFile;
        $oldStat = substr $oldFile, 0, 4;
        `mv $newFile $oldFile`;
        $line = `grep -H $newFile ../rinexantls/*`;
        @_ = split(':', $line);
	$rnxantls = $_[0];
        `mv $rnxantls temp`;
	open NEWRAS, ">$rnxantls";
        open TMP, 'temp';
        while (<TMP>) {
            if (/^$newFile/) {
                s/$newFile/$oldFile/;
            }
            print NEWRAS $_;
        }
        close NEWRAS;
        close TMP;
    }
}
unlink 'temp';
