# rf2m

**rf2m** is a software developed to modify genomes based on VCF files with SNPs and indels present in cell lines or strains. The software also creates a new annotation file (GTF) for the mutant genome with updated coordinates of the genes.

## Requirements

To run **rf2m** you only need to have Perl (I think any versions works, but my Perl version is 5.26.1) installed on your machine. The software also does not need to be installed.

## FASTA file

The rf2m was designed to process fasta files from [ENSEMBL database](https://www.ensembl.org/index.html).

## VCF file format

The specifications for VCF format are described [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf). The script does not work for structural variants, only for SNPs and small indels (there isn't a maximum number of nucleotides in a INDEL supported by our tool. The software will process as many nucleotides as the REF/ALT field have. However, REF/ALT field must have only letters (A,C,G,T,a,c,g,t)). In the present version of **rf2m** only the SNPs and indels with **PASS** in the FILTER field will be processed. We are studying to include an option where the user can choose the Filter value in future releases. If the nucleotides in REF column in VCF file does not match with the nucleotides in the same positions on genome, the script will skip this mutation. When there are overlapping indels only the first will be processed. The script also does not process the lines in VCF with more than one alternative allele in ALT field. Below has an example: The first line will be processed, and the second will not be processed.

|CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO |
|------|-----|----|-----|-----|------|--------|------|
|20    |2    |.   |TC   |TCA  |.     |PASS    |DP=100|
|20    |2    |.   |TC   |TG,T |.     |PASS    |DP=100|

## Running

Run the command lines below in the folder where you downloaded **rf2m**:

`chmod +x format_vcf.pl`

`chmod +x genome_creator.pl`

`chmod +x gtf_creator.pl`

`./format_vcf.pl --vcf input.vcf`

`./genome_creator.pl --genome input_genome.fasta --outGenome output_genome.fasta`

`./gtf_creator.pl --gtf input.gtf --outGTF output.gtf`

### Options
|Option      |Description                                                                       |
|------------|----------------------------------------------------------------------------------|
|--vcf       |File with SNPs and indels in VCF format                                           |
|--genome    |File with reference genome in FASTA format                                        |
|--gtf       |File with annotation for reference genome in GTF format                           |
|--outGenome |Name for the output file of the modified genome in FASTA format                   |
|--outGTF    |Name for the output file of the annotation for the modified genome in GTF format  |
