"""
Date: April 11, 2022

Author: Karine Choquet

This script will analyze processing of human mitochondrial transcripts.

Usage: python mito_RNA_processing_analysis_for_github.py $gtf $bam $out_3prime $out_5prime $out_reads

gtf: gtf file containing the genomic coordinates of all mitochondrial transcripts
bam: alignment file including multi-mapping reads
out_3prime: number of processed and unprocessed reads at the 3'-end for each transcript
out_5prime: number of processed and unprocessed reads at the 5'-end for each transcript
out_reads: 3' and 5' processing status for each read

"""




import sys
import numpy as np
import pandas as pd
import pysam
from collections import Counter

import re
import math

import pybedtools
from pybedtools import BedTool


# CONFIG

hela_gtf = pd.read_table(sys.argv[1], header=None, sep="\t", dtype={'chrom':str, 'build':str, 'feature':str, 'start':int, 'end':int, 'score1':int, 'strand':str, 'score2':int, 'info':str})
#hela_gtf = pd.read_table(sys.argv[1], header=None, sep="\t")
hela_gtf.columns = ['chrom','build','feature','start','end','score1','strand','score2','info']
hela_gtf_mito = hela_gtf[hela_gtf['chrom']=='h_MT'].reset_index(drop=True)


bamFile = pybedtools.BedTool(sys.argv[2])


# Parameters
polyA_window = 15
TSS_window_post = 50
TSS_window_pre = 15

min_overlap = 15
processing_window = 15


# Functions

# Get read ends and read starts

def get_mito_features_bedtool_for_read_ends_starts(annotation_df, polyA_window, TSS_window_pre, TSS_window_post):

    # make a set for all 3'SS coordinates
    features = []

    # loop through a file with intron coordinates
    # check if feature is an exon
    for i in range(0,len(annotation_df)):
        feature = annotation_df['feature'].iloc[i]   # feature
        chrom = annotation_df['chrom'].iloc[i] # chromosome
 
        if (feature == 'transcript') and (chrom == 'h_MT'):
            start = int(annotation_df['start'].iloc[i])             # start coordinate of intron (last base of exon)
            end = int(annotation_df['end'].iloc[i])                 # end coordinate of intron (last base of intron)
            gene = annotation_df['info'].iloc[i].split(";")[4].split('"')[1] # gene name
            gene_biotype = annotation_df['info'].iloc[i].split(";")[6].split('"')[1]
            strand = annotation_df['strand'].iloc[i]                # strand of gene with intron
   
            
            if gene_biotype in ['protein_coding','Mt_rRNA','Mt_tRNA']:
                if (strand=='+'):
                    TSS_start = start - TSS_window_pre
                    TSS_end = start + TSS_window_post
                    
                    polyA_start = end - polyA_window
                    polyA_end = end + polyA_window
                
                if (strand=='-'):
                    TSS_start = end - TSS_window_post
                    TSS_end = end + TSS_window_pre
                    
                    polyA_start = start - polyA_window
                    polyA_end = start + polyA_window

                features.append([chrom,str(polyA_start),str(polyA_end),gene,'3prime_end',strand,gene_biotype])
                features.append([chrom,str(TSS_start),str(TSS_end),gene,'5prime_end',strand,gene_biotype])
                features.append([chrom,start,end,gene,'gene_body',strand,gene_biotype])
                
            elif gene_biotype == "origin":
                features.append([chrom,start,end,gene,'origin',strand,gene_biotype])
              
    features.append(['h_MT',1,16569,'h_MT_+','intergenic','+','intergenic'])
    features.append(['h_MT',1,16569,'h_MT_-','intergenic','-','intergenic'])
    features_bedtool = BedTool(features)
    return features_bedtool



def get_read_end_bedtool(bamFile):

    bedFile = bamFile.bam_to_bed()
    bedFile_df = bedFile.to_dataframe()
        
    read_end = []
        
    for i in range(0,len(bedFile_df)):

        chrom = str(bedFile_df['chrom'].iloc[i])
        start = bedFile_df['start'].iloc[i]
        end = bedFile_df['end'].iloc[i]
        read = bedFile_df['name'].iloc[i]
        score = bedFile_df['score'].iloc[i]
        strand = bedFile_df['strand'].iloc[i]

        if (strand == "-"):
            pos_1 = start
            pos_2 = start + 1

        if (strand == "+"):
            pos_1 = end - 1
            pos_2 = end

        read_end.append([chrom,str(pos_1),str(pos_2),read,str(score),strand])

    read_end_bedtool = BedTool(read_end)
    return read_end_bedtool



def get_read_start_bedtool(bamFile):

    bedFile = bamFile.bam_to_bed()
    bedFile_df = bedFile.to_dataframe()
        
    read_start = []
        
    for i in range(0,len(bedFile_df)):

        chrom = str(bedFile_df['chrom'].iloc[i])
        start = bedFile_df['start'].iloc[i]
        end = bedFile_df['end'].iloc[i]
        read = bedFile_df['name'].iloc[i]
        score = bedFile_df['score'].iloc[i]
        strand = bedFile_df['strand'].iloc[i]

        if (strand == "-"):
            pos_1 = end - 1
            pos_2 = end

        if (strand == "+"):
            pos_1 = start
            pos_2 = start + 1

        read_start.append([chrom,str(pos_1),str(pos_2),read,str(score),strand])

    read_start_bedtool = BedTool(read_start)
    return read_start_bedtool




def get_intersect_read_ends_starts(read_ends, intron_info):

    intersect = read_ends.intersect(intron_info, wo=True, s=True) # intersect reads from bam file with 3' splice site coordinates, ensure strandedness
    intersect_df = intersect.to_dataframe(names=['chrom_read', 'start_read', 'end_read', 'name_read', 'qual_read', \
                                           'strand_read', 'chr_feature', 'start_feature', \
                                           'end_feature', 'name_gene', 'name_feature', 'strand_feature', 'biotype_feature', 'count'], \
                               dtype={"chrom_read": str, "start_read": int, "end_read": int, \
                                     "name_read": str, "qual_read": int, "strand_read": str, \
                                    "chr_feature": str, "start_feature": int, "end_feature": int, "name_gene": str, \
                                     "name_feature": str,"strand_feature": str, "biotype_feature":str, "count": int}) # convert to a dataframe

    return intersect_df



def get_read_end_mapping(intersect_df):

    read_ends = {}

    for i in range(0,len(intersect_df)):

        # get read name and feature type
        read = intersect_df['name_read'].iloc[i]
        feature = intersect_df['name_feature'].iloc[i]
        biotype = intersect_df['biotype_feature'].iloc[i]
        strand = intersect_df['strand_read'].iloc[i]
        
        feature_type = feature + "__" + biotype + "," + strand

        # check if read name is in the dictionary, if not save it
        if read not in read_ends.keys():

            # make a new dictionary for the read and end mapping info
            read_ends[read] = [feature_type]

        # check if read name is in the dictionary, if not save it
        if read in read_ends.keys():

            # if end mapping info is different, append it to the dictionary
            if (feature_type not in read_ends[read]):
                read_ends[read].append(feature_type)
    
    return read_ends



def get_read_end_stats(read_ends):
    
    read_features = []

    for k, v in read_ends.items():

        if (len(v) == 1):
            read_features.append([k,v[0]])

        if (len(v) > 1):
            if ("3prime_end__Mt_rRNA,+" in v):
                read_features.append([k,"3prime_end__Mt_rRNA,+"])
                
            elif ("3prime_end__Mt_tRNA,+" in v):
                read_features.append([k,"3prime_end__Mt_tRNA,+"])
                
            elif ("3prime_end__protein_coding,+" in v):
                read_features.append([k,"3prime_end__protein_coding,+"])
                
            elif ("3prime_end__Mt_rRNA,-" in v):
                read_features.append([k,"3prime_end__Mt_rRNA,-"])
                
            elif ("3prime_end__Mt_tRNA,-" in v):
                read_features.append([k,"3prime_end__Mt_tRNA,-"])
                
            elif ("3prime_end__protein_coding,-" in v):
                read_features.append([k,"3prime_end__protein_coding,-"])
                
            elif ("gene_body__Mt_rRNA,+" in v):
                read_features.append([k,"gene_body__Mt_rRNA,+"])
                
            elif ("gene_body__Mt_tRNA,+" in v):
                read_features.append([k,"gene_body__Mt_tRNA,+"])
                
            elif ("gene_body__protein_coding,+" in v):
                read_features.append([k,"gene_body__protein_coding,+"])
                
            elif ("gene_body__Mt_rRNA,-" in v):
                read_features.append([k,"gene_body__Mt_rRNA,-"])
                
            elif ("gene_body__Mt_tRNA,-" in v):
                read_features.append([k,"gene_body__Mt_tRNA,-"])
                
            elif ("gene_body__protein_coding,-" in v):
                read_features.append([k,"gene_body__protein_coding,-"])
            else:
                read_features.append([k,"undetermined"])

    read_features_df = pd.DataFrame(read_features)
    read_features_df.columns = ['read','end_feature']
    
    return read_features_df



def get_read_start_stats(read_ends):
    
    read_features = []

    for k, v in read_ends.items():

        if (len(v) == 1):
            read_features.append([k,v[0]])

        if (len(v) > 1):
            if ("5prime_end__protein_coding,+" in v):
                read_features.append([k,"5prime_end__protein_coding,+"])
                
            elif ("5prime_end__protein_coding,-" in v):
                read_features.append([k,"5prime_end__protein_coding,-"])
                
            elif ("gene_body__protein_coding,+" in v):
                read_features.append([k,"gene_body__protein_coding,+"])
                
            elif ("gene_body__protein_coding,-" in v):
                read_features.append([k,"gene_body__protein_coding,-"])
            else:
                read_features.append([k,"undetermined"])

    read_features_df = pd.DataFrame(read_features)
    read_features_df.columns = ['read','start_feature']
    
    return read_features_df

###

hela_mito_bedtool = get_mito_features_bedtool_for_read_ends_starts(hela_gtf_mito, polyA_window, TSS_window_pre, TSS_window_post)

# get read ends and turn into a bedtool for intersecting 
read_ends = get_read_end_bedtool(bamFile)

# intersect read ends with genome features
intersect_ends = get_intersect_read_ends_starts(read_ends, hela_mito_bedtool)

# get read ends dictionary
read_end_mapping = get_read_end_mapping(intersect_ends)

# get read end mapping statistics
read_end_stats = get_read_end_stats(read_end_mapping)


# get read starts and turn into a bedtool for intersecting 
read_starts = get_read_start_bedtool(bamFile)

# intersect read starts with genome features
intersect_starts = get_intersect_read_ends_starts(read_starts, hela_mito_bedtool)

# get read starts dictionary
read_end_mapping_starts = get_read_end_mapping(intersect_starts)

# get read end mapping statistics
read_start_stats = get_read_start_stats(read_end_mapping_starts)




#####

## Processing intermediates

def get_mito_features_bedtool_cov(annotation_df, polyA_window, processing_window):

    # make a set for all 3'SS coordinates
    features = []

    # loop through a file with intron coordinates
    # check if feature is an exon
    for i in range(0,len(annotation_df)):
        feature = annotation_df['feature'].iloc[i]   # feature
        chrom = annotation_df['chrom'].iloc[i] # chromosome
 
        if (feature == 'transcript') and (chrom == 'h_MT'):
            start = int(annotation_df['start'].iloc[i])             # start coordinate of intron (last base of exon)
            end = int(annotation_df['end'].iloc[i])                 # end coordinate of intron (last base of intron)
            gene = annotation_df['info'].iloc[i].split(";")[4].split('"')[1] # gene name
            gene_biotype = annotation_df['info'].iloc[i].split(";")[6].split('"')[1]
            strand = annotation_df['strand'].iloc[i]                # strand of gene with intron
   
            
            if gene_biotype in ['protein_coding','Mt_rRNA','Mt_tRNA']:
                if (strand=='+'):
                    polyA_start = end - polyA_window
                    polyA_end = end + polyA_window
                    
                    prefix_start = start - (processing_window * 4)
                    prefix_end = start - processing_window
                    
                    suffix_start = end + processing_window
                    suffix_end = end + (processing_window * 4)
                
                if (strand=='-'):
                    polyA_start = start - polyA_window
                    polyA_end = start + polyA_window
                    
                    prefix_start = end + processing_window
                    prefix_end = end + (processing_window * 4)
                    
                    suffix_start = start - (processing_window * 4)
                    suffix_end = start - processing_window

                features.append([chrom,str(polyA_start),str(polyA_end),gene,'3prime_end,'+gene_biotype,strand])
                features.append([chrom,str(prefix_start),str(prefix_end),gene,'preGeneBody,'+gene_biotype,strand])
                features.append([chrom,str(suffix_start),str(suffix_end),gene,'postGeneBody,'+gene_biotype,strand])
                features.append([chrom,start,end,gene,gene_biotype,strand])
                
            elif gene_biotype == "origin":
                features.append([chrom,start,end,gene,gene_biotype,strand])
              
    features.append(['h_MT',1,16569,'h_MT_+','intergenic','+'])
    features.append(['h_MT',1,16569,'h_MT_-','intergenic','-'])
    features_bedtool = BedTool(features)
    return features_bedtool



# function to create a dataframe with reads that span transcripts
def get_transcript_intersect(transcript_df, bam_file):
    # get reads that span 3' splice sites and convert to a dataframe
    bedFile = bam_file.bam_to_bed(cigar=True, tag='NM') # convert bam file to bed file, keep cigar string and NM (edit distance) tag
    intersect = bedFile.intersect(transcript_df, wo=True, s=True) # intersect reads from bam file with 3' splice site coordinates, ensure strandedness
    df = intersect.to_dataframe(names=['chr_aln', 'start_aln', 'end_aln', 'name_aln', 'qual_aln', \
                                           'strand_aln', 'cigar_aln', 'chr_transcript', 'start_transcript', \
                                           'end_transcript', 'name_gene', 'biotype_gene', 'strand_gene', 'count'], \
                               dtype={"chr_aln": str, "start_aln": int, "end_aln": int, \
                                     "name_aln": str, "qual_aln": int, "strand_aln": str, \
                                     "cigar_aln": str, "chr_transcript": str, "start_transcript": int, \
                                     "end_transcript": int, "name_gene": str, \
                                     "biotype_gene": str,"strand_gene": str, "count": int}) # convert to a dataframe
    return df



# every read that spans a transcript in the dataset
def get_read_junctions_dictionary(intersect_df, min_overlap):

    read_junctions = {}

    

    for i in range(0,intersect_df.shape[0]):       

        # define the read name
        read_name = intersect_df['name_aln'].iloc[i]
        gene_name = intersect_df['name_gene'].iloc[i]
        chrom = intersect_df['chr_transcript'].iloc[i]
        transcript_start = intersect_df['start_transcript'].iloc[i]
        transcript_end = intersect_df['end_transcript'].iloc[i]
        biotype = intersect_df['biotype_gene'].iloc[i]
        strand = intersect_df['strand_gene'].iloc[i]
        read_overlap = intersect_df['count'].iloc[i]
        
        # set variables for parsing the cigar string
        pattern = re.compile('([MIDNSHPX=])')
        Consumes_Query = ["M", "I", "S", "=", "X"]
        Consumes_Reference = ["M", "D", "N", "=", "X"] 
        
        # parse cigar string into a list of tuples for easy parsing
        Sep_Values = pattern.split(intersect_df.cigar_aln[i])[:-1]
        CigarPairs = list((Sep_Values[n:n+2] for n in range(0, len(Sep_Values), 2)))
        
        # get the 3' softclip length
        if intersect_df.strand_aln[i]=="+":        
            last=len(CigarPairs)
            if(CigarPairs[last-1][1]=='S'):
                clip_3prime=CigarPairs[last-1][0]
            elif(CigarPairs[last-1][1]=='H'):
                clip_3prime=CigarPairs[last-1][0]
            elif(CigarPairs[last-1][1]!='S' or CigarPairs[last-1][1]!='H'):
                clip_3prime=0
                
        if intersect_df.strand_aln[i]=="-":
            if(CigarPairs[0][1]=='S'):
                clip_3prime=CigarPairs[0][0]
            elif(CigarPairs[0][1]=='H'):
                clip_3prime=CigarPairs[0][0]
            elif(CigarPairs[0][1]!='S' or CigarPairs[0][1]!='H'):
                clip_3prime=0
                
        # get the 5' softclip length
        if intersect_df.strand_aln[i]=="+":        
            if(CigarPairs[0][1]=='S'):
                clip_5prime=CigarPairs[0][0]
            elif(CigarPairs[0][1]=='H'):
                clip_5prime=CigarPairs[0][0]
            elif(CigarPairs[0][1]!='S' or CigarPairs[0][1]!='H'):
                clip_5prime=0
                
        if intersect_df.strand_aln[i]=="-":
            last=len(CigarPairs)
            if(CigarPairs[last-1][1]=='S'):
                clip_5prime=CigarPairs[last-1][0]
            elif(CigarPairs[last-1][1]=='H'):
                clip_5prime=CigarPairs[last-1][0]
            elif(CigarPairs[last-1][1]!='S' or CigarPairs[last-1][1]!='H'):
                clip_5prime=0
        
        # set up variables for measuring the length of cigar string operators
        CigarOp_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        start_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        end_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        intron_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        currentloc = int(intersect_df['start_aln'].iloc[i])
        
        
        # go through list of cigar strings and grab splicing information
        for cigar_Entry in CigarPairs:

            op_Length = int(cigar_Entry[0]) # get length of cigar operator
            cigarOp = cigar_Entry[1] # get type of cigar operator  
            CigarOp_counts[cigarOp] += op_Length # add the cigar operator length to the counts dictionary
            cigarOp_start=currentloc # get the starting coordinate of the cigar operator

            if (cigarOp in Consumes_Reference):
                currentloc=currentloc+op_Length # add the cigar operator length to the current location coordinate 

            cigarOp_end=currentloc # get the ending coordinate of the cigar operator

            # gather information if the portion of the cigar string spans the designated intron start
            if (cigarOp_start<transcript_start-min_overlap and cigarOp_end>=transcript_start-min_overlap):
                if (cigarOp_end>=transcript_start+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(transcript_start-min_overlap)+1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=transcript_start-min_overlap and cigarOp_end<transcript_start+min_overlap):
                count=op_Length
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary       

            elif (cigarOp_start<transcript_start+min_overlap and cigarOp_end>=transcript_start+min_overlap):
                if (cigarOp_start<=transcript_start-min_overlap):
                    count=min_overlap*2
                else:
                    count=(transcript_start+min_overlap)-cigarOp_start-1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string is within the intron
            if (cigarOp_start<transcript_start and cigarOp_end>=transcript_start):
                if (cigarOp_end>=transcript_end):
                    count=transcript_end-transcript_start
                else:
                    count=cigarOp_end-transcript_start
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=transcript_start and cigarOp_end<transcript_end):
                count=op_Length
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<transcript_end and cigarOp_end>=transcript_end):
                if (cigarOp_start<=transcript_start):
                    count=transcript_end-transcript_start
                else:
                    count=transcript_end-cigarOp_start
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string spans the designated intron end
            if (cigarOp_start<transcript_end-min_overlap and cigarOp_end>=transcript_end-min_overlap):
                if (cigarOp_end>=transcript_end+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(transcript_end-min_overlap)
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=transcript_end-min_overlap and cigarOp_end<transcript_end+min_overlap):
                count=op_Length
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<transcript_end+min_overlap and cigarOp_end>=transcript_end+min_overlap):
                if (cigarOp_start<=transcript_end-min_overlap):
                    count=min_overlap*2
                else:
                    count=(transcript_end+min_overlap)-cigarOp_start
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary
        
        
        # Remove reads that have long portions that are "spliced" between two genes (at least for now)
        if (end_counts['N']==0 and start_counts['N']==0):
            if read_overlap > min_overlap:

                # check if read name is in the dictionary, if not save it
                if read_name not in read_junctions.keys():
                    read_junctions[read_name]= [[chrom, transcript_start, transcript_end, gene_name, biotype, strand, read_overlap, clip_5prime, clip_3prime]]

                # check if read name is in the dictionary, if it is proceed to gene information
                elif read_name in read_junctions.keys():
                    read_junctions[read_name].append([chrom, transcript_start, transcript_end, gene_name, biotype, strand, read_overlap, clip_5prime, clip_3prime])


    return read_junctions



def get_read_transcript_mapping_with_processing(read_junctions):
    
    read_property = []
    
    gene_counts = {}

    # loop through all reads in the dictionary
    for read in read_junctions.keys():
        
        for row in read_junctions[read]:

            # characterize the number of spliced and unspliced introns in the read
            gene_status = row[3]
            biotype = row[4]
            strand = row[5]
        
            read_property.append([read, gene_status, biotype, strand])

    read_transcript_df = pd.DataFrame(read_property)
    read_transcript_df.columns = ['read','gene_name', 'biotype', 'strand']

    return read_transcript_df


####

hela_mito_bedtool_cov = get_mito_features_bedtool_cov(hela_gtf_mito, polyA_window, processing_window)

# count reads mapping to whole gene/transcript
# get reads that intersect transcripts
intersect_cov = get_transcript_intersect(hela_mito_bedtool_cov, bamFile)

# get dictionary with all transcripts that a read spans
splice_dictionary = get_read_junctions_dictionary(intersect_cov, min_overlap)

# get gene mapping per read with processing intermediates
process_df = get_read_transcript_mapping_with_processing(splice_dictionary)

# Filter for protein-coding genes
process_df = process_df[process_df['biotype'].str.contains('protein_coding')].drop_duplicates().reset_index(drop=True)

# Add a column indicating presence
process_df['presence'] = 'yes'

# Combine read and gene name
process_df['read,gene_name'] = process_df['read'] + ',' + process_df['gene_name']

# Pivot table
process_piv = process_df.pivot(index='read,gene_name', columns='biotype', values='presence').fillna("no").reset_index()

# Separate read and gene name
process_piv['read'] = process_piv['read,gene_name'].str.split(",").str[0]
process_piv['gene_name'] = process_piv['read,gene_name'].str.split(",").str[1]

# Merge with read end features
process_read_ends = process_piv.merge(read_end_stats, on='read').merge(read_start_stats, on='read')

# Extract the patterns of interest (must cover gene body)
process_read_ends = process_read_ends[(process_read_ends['protein_coding']=='yes')].reset_index(drop=True)

# Add a category based on where the read maps
process_read_ends['category_3prime'] = 'other'

process_read_ends.loc[(process_read_ends['protein_coding']=='yes') & 
                            (process_read_ends['end_feature'].str.startswith('3prime_end__protein')) &
                            (process_read_ends['postGeneBody,protein_coding']=='no'),'category_3prime'] = 'processed_3prime'

process_read_ends.loc[(process_read_ends['protein_coding']=='yes') & 
                            (process_read_ends['postGeneBody,protein_coding']=='yes'),'category_3prime'] = 'unprocessed_3prime'

process_read_ends.loc[(process_read_ends['protein_coding']=='yes') & 
                            (process_read_ends['end_feature'].str.startswith('gene_body')) &
                            (process_read_ends['postGeneBody,protein_coding']=='no'),'category_3prime'] = 'transcribing'


process_read_ends['category_5prime'] = 'other'

process_read_ends.loc[(process_read_ends['protein_coding']=='yes') & 
                            (process_read_ends['start_feature'].str.startswith('5prime_end__protein')) &
                            (process_read_ends['preGeneBody,protein_coding']=='no'),'category_5prime'] = 'processed_5prime'

process_read_ends.loc[(process_read_ends['protein_coding']=='yes') & 
                            (process_read_ends['preGeneBody,protein_coding']=='yes'),'category_5prime'] = 'unprocessed_5prime'

process_read_ends.to_csv(sys.argv[-1], sep="\t", header=True, index=False)


# Count
process_counts_3prime = process_read_ends[(process_read_ends['category_3prime']!='other') & 
                                                   (~process_read_ends['gene_name'].str.contains("anti"))].groupby(['gene_name','category_3prime'])['read'].count().reset_index().pivot(index='gene_name', columns='category_3prime', values='read').fillna(0).reset_index()

process_counts_5prime = process_read_ends[(process_read_ends['category_5prime']!='other') & 
                                                   (~process_read_ends['gene_name'].str.contains("anti"))].groupby(['gene_name','category_5prime'])['read'].count().reset_index().pivot(index='gene_name', columns='category_5prime', values='read').fillna(0).reset_index()


# Write output files
process_counts_3prime.to_csv(sys.argv[3], sep="\t", header=True, index=False)
process_counts_5prime.to_csv(sys.argv[4], sep="\t", header=True, index=False)


