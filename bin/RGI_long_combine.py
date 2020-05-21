#!/usr/bin/env python3

import sys
import csv


def rgi_long_combine(rgi_perf_file, long_file, combined_output):
    # Get the desired information from the RGI file
    with open(rgi_perf_file, 'r') as rgifile:

        # Define a dictionary where desired formatted rgi entries (gene column format) are keys and the items are the number of hits and default gene fraction percentage
        rgi_perf_dict = {}

        # Go through the file and fill out the dictionary using long format
        reader = csv.reader(rgifile, delimiter=',')
        next(reader)
        for row in reader:
            aro_name = "RGI|" + row[2] + "|" + row[0]
            rgi_perf_dict[aro_name] = [row[3], 80]

        # Close the file
        rgifile.close()

    # Get the name of the sample
    sample_name = rgi_perf_file.split("_")[0]

    # Get counts from the provided long_file
    with open(long_file, 'r') as long_file:

        # Define a dictionary where genes are keys and the items are the sample name, number of hits, and gene fraction percentage
        long_dict = {}

        # Go through the file and fill out the dictionary using long format lines
        long_reader = csv.reader(long_file, delimiter=',')
        header = next(long_reader)
        for long_row in long_reader:
            split_gene_name = long_row[1].split("|")
            if split_gene_name[len(split_gene_name)-1] != "RequiresSNPConfirmation":
                long_dict[long_row[1]] = [long_row[0], long_row[2], long_row[3]]

        # Close the file
        long_file.close()


    # Write combined output to the provided output file
    with open(combined_output, 'w', newline='\n') as combined_file:
        combined_write = csv.writer(combined_file, delimiter=',')
        combined_write.writerow(header)
        # Write to the combined file using the dictionaries we created previously
        for rgi_key, rgi_item in rgi_perf_dict.items():
            temp_rgi_write = rgi_item.copy()
            temp_rgi_write.insert(0, rgi_key)
            temp_rgi_write.insert(0, sample_name)
            combined_write.writerow(temp_rgi_write)

        for long_key, long_item in long_dict.items():
            temp_long_write = long_item.copy()
            temp_long_write.insert(1, long_key)
            combined_write.writerow(temp_long_write)


if __name__ == '__main__':
    rgi_long_combine(sys.argv[1], sys.argv[2], sys.argv[3])
