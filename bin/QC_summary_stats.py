import os
import sys
import gzip
import argparse
import glob
import sys
import gzip
import csv
import pandas as pd
import numpy




def parse_cmdline_params(cmdline_params):
	info = "Removes duplicate FASTQ entries from a FASTQ file"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument('-i', '--input_files', nargs='+', required=True,
        	                help='Use globstar to pass a list of sequence files, (Ex: *.fastq.gz)')
	return parser.parse_args(cmdline_params)


def pull_Phred(fastq_files):

    for f in fastq_files: # iterate through each fastq file
        Plist=[]
        Qlist=[]
        Seqlen_list=[]
        num_reads = 0

        fp = open(f, 'r') # open each fastq file;   gzip.open if .gz files
        for line in fp:   # iterate through lines of fastq file
            Ordqual=[]
            Q=[]
            P=[]
            read_id = line
            seq = fp.next()

            #seq = seq[10:len(seq)] # Let's not chop off the umi here since we would be checking the quality after UMI removal and not all samples have UMIs
            
            Seqlen_list.append(len(seq))
            #newseq = seq + spacesep + UMI
            plus = fp.next()
            qual = fp.next()



            for i in range(len(qual)-1): #Exclude the return character
                Ordqual.append(ord(qual[i]))
                Q.append(Ordqual[i]-33)
                P.append(10**(-Q[i]/10))

            Qlist.append(numpy.mean(Q))
            Plist.append(numpy.mean(P))
                
            num_reads += 1
                
        print(f,"mean_probability_nucleotide_error",numpy.mean(Plist))
        print(f,"mean_phred_score",numpy.mean(Qlist))
        print(f,"total_reads",num_reads)
        print(f,"mean_read_length",numpy.mean(Seqlen_list))

        fp.close()








def print_dict(dict):
  # iterate through UMIs and repeat counts and print those
  dups= "Repeat UMI Error"
  for k, v in dict.items():
	  if v != 1:
		  print k, dups, v
	  else:
	  	print k, v

def print_Rep(dict):
  # iterate through UMIs and repeat counts and print those
  dups= "Repeat UMI Error"
  for k, v in dict.items():
	  if len(v) != 1:
		  print k, dups, v
	  else:
	  	print k, v


#Apply the previous functions; Print UMIs and counts
if __name__ == "__main__":
  opts = parse_cmdline_params(sys.argv[1:])
  fastq_files = opts.input_files
  #read_dic = pull_UMI(fastq_files)
  #print_dict(read_dic[0])
  #print_Rep(read_dic[1])
  #print ReWriter(read_dic[1])
  pull_Phred(fastq_files)
