#!/usr/bin/perl -w

#Script to calculate the GC content of genomic sequences. Adopted script (get_gc_content.pl) from Jennifer Meneghin .
####################################################################################################
### Get GC Content                                                                               ###
#Usage : perl gc-calculator.pl <input-fasta-file> <output-file-path>                             ###
####################################################################################################
#---------------------------------------------------------------------------------------------------------------------------
#Deal with passed parameters
#---------------------------------------------------------------------------------------------------------------------------
if ($#ARGV == -1) {
    usage();
    exit;
}
$fasta_file = $ARGV[0];
$out_file = $ARGV[1];
unless ( open(IN, "$fasta_file") ) {    
    print "Got a bad fasta file: $fasta_file\n\n";
    exit;
}
unless ( open(OUT, ">$out_file") ) {
    print "Couldn't create $out_file\n";
    exit;
}
print "Parameters:\nfasta file = $fasta_file\noutput file = $out_file\n\n";
#---------------------------------------------------------------------------------------------------------------------------
#The main event
#---------------------------------------------------------------------------------------------------------------------------
#print OUT "ID\t% GCContent\tTotal Count\tG Count\tC Count\tA Count\tT Count\n";
$seq = "";
while (<IN>) {
    chomp;
    if (/^>/) {
	#finish up previous line.
	if (length($seq) > 0) {
	    &process_it;
	}
	#start new line.
	$id = $_;
	$id =~ s/^>(.+?)\s.+$/$1/g;
	print OUT "$id\t";
    }
    else {
	$seq = $seq . $_;
    }
}

#finish up last line.
&process_it;

close(IN);
close(OUT);

sub usage {
    print "Get GC Content\n";
    print "Usage: gc-calculator.pl <input-fasta-file> <output-file-path>  \n";
    print "This program takes a fasta file as it's first parameter and output file path as second.\n\n";
    print "It returns a tab delimited file (output-file-path): column 1 = header ID (everything between \">\"\n";
    print "and the first space in the header), and column 2 = gc content for the fasta entry.\n\n";
    print "Adopted from script by Jennifer Meneghin\n";
    print "July 23, 2009\n\n";

}

sub process_it {
    @letters = split(//, $seq);
    $gccount = 0;
    $totalcount = 0;
    $acount = 0;
    $tcount = 0;
    $gcount = 0;
    $ccount = 0;
    foreach $i (@letters) {
	if (lc($i) =~ /[a-z]/) {
	    $totalcount++;
	}
	if (lc($i) eq "g" || lc($i) eq "c") {
	    $gccount++;
	}
	if (lc($i) eq "a") {
	    $acount++;
	}
	if (lc($i) eq "t") {
	    $tcount++;
	}
	if (lc($i) eq "g") {
	    $gcount++;
	}
	if (lc($i) eq "c") {
	    $ccount++;
	}
    }
    if ($totalcount > 0) {
	$gccontent = (100 * $gccount) / $totalcount;
    }
    else {
	$gccontent = 0;
    }
    print OUT "$gccontent\n";
    $seq = "";
}
