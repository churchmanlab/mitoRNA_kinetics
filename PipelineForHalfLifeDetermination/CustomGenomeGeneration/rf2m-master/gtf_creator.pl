#!/usr/bin/perl

use Getopt::Long;

$gtf = '';
$output_gtf = '';

GetOptions('gtf=s' => \$gtf, 'outGTF=s' => \$output_gtf) or die "Usage: $0 --gtf NAME --outGTF NAME\n";

$time = localtime;
print "Initiating GTF creation step [$time]\n";

%index;
open(GTF, "$gtf") or die "Failed to open $gtf\n";
while ($line_gtf = <GTF>) {
    chomp($line_gtf);    	
    if ($line_gtf =~ m/^\#/) {
	next;
    }
    @gtf_tab = split("\t", $line_gtf);
    $index{$gtf_tab[0]}{$gtf_tab[3]} = $gtf_tab[3];
    $index{$gtf_tab[0]}{$gtf_tab[4]} = $gtf_tab[4];
}
close(GTF);

undef $line_gtf;
undef @gtf_tab;

open(VCF, "mutations_keeped.vcf") or die "Cannot open file mutations_keeped.vcf. This file does not exist. Run genome_creator.pl first!\n";
while ($vcf_line = <VCF>) {
    if ($vcf_line =~ m/^\#/) {
	next;
    }
    @vcf_tab = split("\t",$vcf_line);
    $chr = $vcf_tab[0];
    $pos = $vcf_tab[1];
    $ref = $vcf_tab[3];	    
    $alt = $vcf_tab[4];
    $diff = length($alt) - length($ref);
    if ($diff != 0) {
	foreach $position (keys %{$index{$chr}}) {
	    if ($index{$chr}{$position} >= $pos) {
		$index{$chr}{$position} += $diff;
	    }
	}
    }
}
close(VCF);

undef $vcf_line;
undef $vcf_tab;

$time = localtime;
print "Writting new coordinates of genes in $output_gtf [$time]\n";

open(OUT_GTF, ">$output_gtf") or die "Failed to create $output_gtf\n";
open(GTF, "$gtf") or die "Failed to open $gtf\n";
while ($line_gtf = <GTF>) {
    chomp($line_gtf);    	
    if ($line_gtf =~ m/^\#/) {
	print OUT_GTF "$line_gtf\n";
	next;
    }
    @gtf_tab = split("\t", $line_gtf);
    $gtf_tab[3] = $index{$gtf_tab[0]}{$gtf_tab[3]};
    $gtf_tab[4] = $index{$gtf_tab[0]}{$gtf_tab[4]};
    $new_gtf_line = join("\t", @gtf_tab);
    print OUT_GTF "$new_gtf_line\n";
}
close(GTF);
close(OUT_GTF);

$time = localtime;
print "gtf_creator.pl complete [$time]\n";
