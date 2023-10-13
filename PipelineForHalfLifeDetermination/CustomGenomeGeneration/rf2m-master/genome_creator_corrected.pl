#!/usr/bin/perl

use Getopt::Long;

$genome = '';
$output_genome = '';

GetOptions('genome=s' => \$genome, 'outGenome=s' => \$output_genome) or die "Usage: $0 --genome NAME --outGenome NAME\n";

$time = localtime;

print "Iniciating analysis [$time]\n";

$time = localtime;

print "Opening reference genome file and parsing VCF file [$time]\n";

%vcf_keeped;

open(OUT_GENOME, ">$output_genome") or die "Failed to create $output_genome\n";
open(GENOME, "$genome") or die "Cannot open file $genome. This file does not exist!\n";
$/ = ">";
<GENOME>;
while ($genome_line = <GENOME>) {
    chomp($genome_line);
    ($header, @sequence) = split("\n",$genome_line);
    $header =~ m/^(.+) 1-/;
    $chr_ref = $1;
    $time = localtime;
    print "Reading chromosome $chr_ref [$time]\n";
    $complete_sequence = join("",@sequence);
    undef @sequence;
    @seq_per_base = split("",$complete_sequence);
    undef $complete_sequence;
    open(VCF, "mutations_accepted.vcf") or die "Cannot open mutation_accepted.vcf. This file does not exist. Run format_vcf.pl first!\n";
    $/ = "\n";
    while ($vcf_line = <VCF>) {
	chomp($vcf_line);
	if ($vcf_line =~ m/^\#/) {
	    next;
	}
	@vcf_tab = split("\t",$vcf_line);
	$chr = $vcf_tab[0];
	$pos = $vcf_tab[1] - 1;
	$ref = uc $vcf_tab[3];	    
	$alt = uc $vcf_tab[4];
	$filter = $vcf_tab[6];
	if ($chr eq $chr_ref) {
	    if (length($ref) > 1) {
		@chg = split("", $ref);
		$i = 0;
		$len = scalar(@chg) - 1;
		for ($i; $i <= $len; $i++) {
		    if ($i == 0) {
			if ($chg[0] eq $seq_per_base[$pos]) {
			    if ($alt eq ".") {
				$seq_per_base[$pos] = "x";
				$changes += 1;
				$pos += 1;
				$vcf_keeped{$vcf_line} = 1;
			    }
			    else {
				$seq_per_base[$pos] = $alt;
				$changes += 1;
				$pos += 1;
				$vcf_keeped{$vcf_line} = 1;
			    }
			}
			else {
			    last;
			}
		    }
		    else {
			$seq_per_base[$pos] = "x";
			$pos += 1;
		    }
		}
		$pos = $vcf_tab[1] - 1;
	    }
	    else {
		if ($ref eq $seq_per_base[$pos]) {
		    if ($alt eq ".") {
			$seq_per_base[$pos] = "x";
			$changes += 1;
			$vcf_keeped{$vcf_line} = 1;
			next;
		    }
		    else {
			$seq_per_base[$pos] = $alt;
			$changes += 1;
			$vcf_keeped{$vcf_line} = 1;
			next;
		    }
		}
		else {
		    next;
		}
	    }
	}
    }
    close(VCF);
    $time = localtime;
    print "Changing reference genome sequences and writing into $output_genome [$time]\n";
    $new_seqs = join("",@seq_per_base);
    undef @seq_per_base;
    $new_seqs =~ s/x//g;
    @new_genome = split("",$new_seqs);
    undef $new_seqs;
    print OUT_GENOME ">$header\n";
    $x = 0;
    $j = 0;
    $len = scalar(@new_genome) - 1;
    for ($x; $x <= $len; $x++) {
	if ($j == 59) {
	    print OUT_GENOME "$new_genome[$x]\n";
	    $j = 0;
	}
	else {
	    if ($x == $len) {
		print OUT_GENOME "$new_genome[$x]\n";
	    }
	    else {
		print OUT_GENOME "$new_genome[$x]";
		$j += 1;
	    }
	}
    }
    undef @new_genome;
    $/ = ">";
}
close(GENOME);

$time = localtime;
print "Writing mutations that were found in reference genome to the file mutations_keeped.vcf [$time]\n";

open(OUT_VCF, ">mutations_keeped.vcf") or die "Failed to create mutations_keeped.vcf\n";
foreach $key (sort(keys(%vcf_keeped))) {
    print OUT_VCF "$key\n";
}
close(OUT_VCF);

print "$changes changes were made\n";
$time = localtime;
print "Analysis complete [$time]\n";
print "Please, proceed to gtf_creator.pl\n";		
