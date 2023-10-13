#!/usr/bin/perl

use Getopt::Long;

$vcf = '';

GetOptions('vcf=s' => \$vcf) or die "Usage $0 --vcf NAME\n";

$time = localtime;

print "Formating $vcf file to genome_creator.pl [$time]\n";

$ini_last = 0;
$fin_last = 0;
$chr = "";

open(VCF, "$vcf") or die "Failed to open $vcf\n";
open(OUT_VCF, ">mutations_accepted.vcf") or die "Failed to create the formatted VCF file\n";
while ($vcf_line = <VCF>) {
    chomp($vcf_line);
    if ($vcf_line =~ m/\#/) {
	print OUT_VCF "$vcf_line\n";
	next;
    }
    else {
	@vcf_tab = split("\t",$vcf_line);
	if ($vcf_tab[6] ne "PASS") {
	    next;
	}
	else {
	    if ($vcf_tab[4] =~ m/\,/) {
		next;
	    }
	    else {
		if (not $vcf_tab[4] =~ m/\A([ACGT]+|\.)\z/i) {
		    next;
		}
		else {
		    $vcf_chr = $vcf_tab[0];
		    if ($vcf_chr eq $chr) {
			$ini = $vcf_tab[1];
			if ($ini >= $ini_last and $ini <= $fin_last) {
			    next;
			}
			else {
			    if ($alt eq ".") {
				$fin = $ini + (length($vcf_tab[3]) - 1);
				$ini_last = $ini;
				$fin_last = $fin;
				print OUT_VCF "$vcf_line\n";
			    }
			    else {
				$fin = $ini + abs(length($vcf_tab[3]) - length($vcf_tab[4]));
				$ini_last = $ini;
				$fin_last = $fin;
				print OUT_VCF "$vcf_line\n";
			    }
			}
		    }
		    else {
			$ini_last = 0;
			$fin_last = 0;
			$chr = $vcf_chr;
			$ini = $vcf_tab[1];
			if ($ini >= $ini_last and $ini <= $fin_last) {
			    next;
			}
			else {
			    if ($alt eq ".") {
				$fin = $ini + (length($vcf_tab[3]) - 1);
				$ini_last = $ini;
				$fin_last = $fin;
				print OUT_VCF "$vcf_line\n";
			    }
			    else {
				$fin = $ini + abs(length($vcf_tab[3]) - length($vcf_tab[4]));
				$ini_last = $ini;
				$fin_last = $fin;
				print OUT_VCF "$vcf_line\n";
			    }
			}
		    }
		}
	    }
	}
    }
}

$time = localtime;

print "Success! Please, proceed to genome_creator.pl [$time]\n"
