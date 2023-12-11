#!/usr/bin/perl
use warnings FATAL => 'all';;
use strict;
use Getopt::Long;

##############################################################
#  script: softmask2hardmask.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.01.20
#  description: identify base converted reads 
#  example: perl mpileup2vcf.pl -v input.vcf -m input.mpleup -s 2 -q 30 -o output.reads.list.txt
##############################################################

my ($vcf, $mpileup, $output, $score_cutoff, $qual_cutoff);
$score_cutoff = 2;
$qual_cutoff = 30;
GetOptions('v|vcf:s' => \$vcf,
	   'mpileup|mpileup:s' => \$mpileup,
	   'o|output:s' => \$output,
	   's|score_cutoff:i' => \$score_cutoff,
	   'q|qual_cutoff:i' => \$qual_cutoff);

my $vcf_fh = read_file($vcf);
my %vcf = parse_vcf_file($vcf_fh);
my %reads = ();
my $mpileup_fh = read_file($mpileup);

while (<$mpileup_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $pos, $ref_allele, $depth, $reads_bases, $reads_quals, $reads_ids) = split /\t/, $_;
    if (exists $vcf{$chr}{$pos}) {
	# print "raw mpileup: $_\n";                                                                                             # remove the start and end marks of reads mapping
	$reads_bases =~ s/\^\S//gi;
	$reads_bases =~ s/\$//gi;
	# print "reads_bases: $reads_bases\n";
	my %indels = ();
	my $j = 0; # indel index         
	while ($reads_bases =~ /((?:\S)(?:\+|\-)[0-9]+[ACGTNacgtn]+)/g) {
	    $j++;
	    my $indel_raw = $&;
	    # my $match_start = $-[0] + 1;
	    # my $match_end = $+[0];
	    $indels{$j}{'indel_raw'} = $&;
	    ($indels{$j}{'lead_base'}, $indels{$j}{'type'}, $indels{$j}{'length'}, $indels{$j}{'indel_adjusted'}) = ($indels{$j}{'indel_raw'} =~ /(\S)(\+|\-)([0-9]+)([ACGTNacgtn]+)/);
	    $indels{$j}{'indel_adjusted'} = substr($indels{$j}{'indel_adjusted'}, 0, $indels{$j}{'length'});
	    if ($indels{$j}{'lead_base'} =~ /(\.|\,)/) {
		$indels{$j}{'lead_base'} = $ref_allele;
	    } else {
		# print "corner case! lead_base = $indels{$j}{'lead_base'}\n";
		$indels{$j}{'lead_base'} = uc $indels{$j}{'lead_base'};
	    }
	    if ($indels{$j}{'type'} eq '+') {
		$indels{$j}{'indel_adjusted'} = $indels{$j}{'lead_base'}. $indels{$j}{'indel_adjusted'};
		$reads_bases =~ s/(?:\S)(?:\+|\-)[0-9]+[ACGTNacgtn]{$indels{$j}{'length'}}/I/;
	    } else {
		$indels{$j}{'indel_adjusted'} = $indels{$j}{'lead_base'};
		$reads_bases =~ s/(?:\S)(?:\+|\-)[0-9]+[ACGTNacgtn]{$indels{$j}{'length'}}/D/;
	    }
	    # print "reads_bases after replacement: $reads_bases\n";
	}
	$j = 0;
	my @reads_ids = split /,/, $reads_ids;
	for (my $i = 1; $i <= $depth; $i++) {
	    my $base = substr $reads_bases, $i - 1, 1;
	    my $qual = substr $reads_quals, $i - 1, 1;
	    my $read_id = $reads_ids[$i-1];
	    $ref_allele = uc $ref_allele;
	    $base = uc $base;
	    $qual = ord($qual) - 33;
	    # if ((($ref_allele eq "T") and ($base eq "C")) or (($ref_allele eq "A") and ($base eq "G"))) {
	    if (($ref_allele eq "T") and ($base eq "C")) {
		if ($qual >= $qual_cutoff) {
		    if (not exists $reads{$read_id}) {
			$reads{$read_id} = 1;
		    } else {
			$reads{$read_id} += 1;
		    }
		}
	    }
	}
    }
}

my $output_fh = write_file($output);

foreach my $id (sort keys %reads) {
    if ($reads{$id} >= $score_cutoff) {
	print $output_fh "$id\n";
    }
}



sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    } else {
        open($fh, $file) or die "can't open $file";
    }
    return $fh;
}

sub write_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  

sub parse_vcf_file {
    my $fh = shift @_;
    my %vcf = ();
    while (<$fh>) {
	chomp;
	 /^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $id, $ref_allele, $alt_allele, $qual, $filter) = split /\t/, $_;
	# if ((($ref_allele eq "T") and ($alt_allele eq "C")) or  (($ref_allele eq "A") and ($alt_allele eq "G"))) {
	$ref_allele = uc $ref_allele;
	$alt_allele = uc $alt_allele;
	if (($ref_allele eq "T") and ($alt_allele eq "C")) {
	    $vcf{$chr}{$pos} = "$ref_allele->$alt_allele";
	}
    }
    return %vcf;
}
	
