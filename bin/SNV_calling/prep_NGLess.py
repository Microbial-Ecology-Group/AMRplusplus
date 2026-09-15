#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 12:35:25 2026

@author: Jared G. Young
"""
from pathlib import Path
import argparse

# version 1.0.0

def main():
    parser = argparse.ArgumentParser(description='Generates a concatenated count matrix and annotation file for the raw files given',
                                     usage='python %(prog)s [-h] [-s save_dir] file_dir fn_struc sample_file')

    parser.add_argument('-aln_wf', type=str, default='Standard', help='Run SNV calling on Standard or Deduplicated aligned reads. Options = Standard (default), Deduped (deduplicated reads).')
    parser.add_argument('-min_len', type=int, default=45, help='Minimum alignment length of reads for SNV calling. Cannot be less than 45 bp.')
    parser.add_argument('-min_ani', type=int, default=97, help='Minimum average nucleotide identity (ANI) of read to reference. Must be ≥ 97% ANI.')
    
    args = parser.parse_args()

    ## Write path inputs and outputs
    if args.aln_wf == "Standard":
        inpath = "test_results/Alignment/SAM_files/Standard/"
        intail = ".alignment.sam"
        outpath = "test_results/Alignment/SNV_BAM_files/SNV_Standard/"
        outtail = ".SNV.unique.sorted.bam"
        out_dir = Path(outpath)
        out_dir.mkdir(parents=True, exist_ok=True) # make outfile path
    elif args.aln_wf == "Deduped":
        inpath = "test_results/Alignment/SAM_files/Deduped/"
        intail = ".alignment.dedup.sam"
        outpath = "test_results/Alignment/SNV_BAM_files/SNV_Deduped/"
        outtail = ".SNV.unique.sorted.dedup.bam"
        out_dir = Path(outpath)
        out_dir.mkdir(parents=True, exist_ok=True) # make outfile path
    else:
        print("Input alignment type" + args.aln_wf + " is not supported")
        
    ## Alignment params 
    if args.min_len < 45:
        print("Minimum alignment length is too short. Must be at least 45bp.")
    minalnlen = args.min_len
    
    if args.min_ani < 97:
        print("Minimum ANI must be ≥ 97%.")
    minani = args.min_ani
    
    ## Constant params 
    unmat = "{unmatch}"   # for unmatch 
    kpif = "{mapped}, {unique}" # for keep if command 
    
    
    ## Write filter.ngl script 
    script_content = f"""
    ngless "1.5"
    import "parallel" version "1.1"
    import "mocat" version "0.0"
    import "samtools" version "0.0"
    
    # A script to filter alignments for SNV calling with metaSNV v2
    
    current = run_for_all(readlines("/test_results/Alignment/SAM_files/sam_filenames.txt"))
    input = samfile('{inpath}' + current + '{intail}')
    
    filtered = select(input) using |mr|:
     mr = mr.filter(min_match_size= {minalnlen} , min_identity_pc= {minani} , action={unmat})
     
    filtered_unique = select(filtered, keep_if=[{kpif}])
    filtered_unique = samtools_sort(filtered_unique)
    write(filtered_unique, ofile='{outpath}' + current + '{outtail}')
    """
    
    # write filter.ngl with user defined params - need to decide where to put this 
    with open("bin/filter.ngl", "w") as f:
        f.write(script_content)
    
    # set path for sam_filenames.txt
    if args.aln_wf == "Standard":
        sam_names = Path("test_results/Alignment/SAM_files/Standard")
    elif args.aln_wf == "Deduped":
        sam_names = Path("test_results/Alignment/SAM_files/Deduped")
        
    # Write filenames.txt
    sam_names_without_ext = [file.stem for file in sam_names.iterdir() if file.is_file()]
        
    with open('test_results/Alignment/SAM_files/sam_filenames.txt','w') as file:
        for item in sam_names_without_ext:
            item = item.split(sep=".")[0]
            file.write(f"{item}\n")
            
    # make file of bam_filepaths.txt
    with open('test_results/Alignment/SAM_files/bam_filepaths.txt','w') as file:
        for item in sam_names_without_ext:
            file.write(f"{outpath}{item}{outtail}\n")

# Define main 
if __name__ == '__main__':
    main()
