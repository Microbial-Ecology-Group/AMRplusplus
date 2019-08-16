#!/usr/bin/env python3

import sys
import csv


def rgi_output(rgi_file):
    # Get the desired information from the RGI file
    with open(rgi_file, 'r') as rgifile:

        # Get the first line of the file and put it into a list
        header = rgifile.readline().split('\t')

        # Get the index of all the elements we want
        category = header.index('Cut_Off')
        best_hit_aro = header.index('Best_Hit_ARO')
        aro = header.index('ARO')
        model_type = header.index('Model_type')

        # Define a dictionary where each best_hit_aro is a key and the items are the rest of the desired elements
        rgi_dict = {}

        # Go through the file and fill out the dictionary without repeats
        reader = csv.reader(rgifile, delimiter="\t")
        for row in reader:
            aro_name = row[best_hit_aro]
            if aro_name not in rgi_dict.keys():
                rgi_dict[aro_name] = [row[aro], row[category], 1, row[model_type]]
            else:
                rgi_dict[aro_name][2] += 1

        # Close the file
        rgifile.close()

    # Get the name of the file to use for the three outputs
    sample_name = rgi_file.split(".")

    perf_sample_name = sample_name
    perf_sample_name.pop()
    perf_sample_name.pop()
    perf_sample_name.insert(1, '_rgi_perfect_hits.csv')
    perf_file_name = ''.join(perf_sample_name)

    strict_sample_name = sample_name
    strict_sample_name.pop()
    strict_sample_name.insert(1, '_rgi_strict_hits.csv')
    strict_file_name = ''.join(strict_sample_name)

    loose_sample_name = sample_name
    loose_sample_name.pop()
    loose_sample_name.insert(1, '_rgi_loose_hits.csv')
    loose_file_name = ''.join(loose_sample_name)

    # Search the dictionary to see which of the three Cut_Offs exist
    perf_in_dict = False
    strict_in_dict = False
    loose_in_dict = False

    for x in rgi_dict.values():
        if x[1] == "Perfect":
            perf_in_dict = True
        if x[1] == "Strict":
            strict_in_dict = True
        if x[1] == "Loose":
            loose_in_dict = True
        # Stop checking if we already know we need to make all three files
        if perf_in_dict and strict_in_dict and loose_in_dict:
            break

    # Write the dictionary to each of the files if cut_off values exist
    # I.e if the file has no perfects, we don't write a file for it
    if perf_in_dict:
        with open(perf_file_name, 'w', newline='\n') as perf_file:
            perf_write = csv.writer(perf_file, delimiter=',')
            perf_write.writerow(["Best_Hit_ARO", "ARO", "Cut_Off",  "Sum_Hits", "Model_Type"])
            for perf_key, perf_item in rgi_dict.items():
                if perf_item[1] == "Perfect":
                    temp_perf_write = perf_item.copy()
                    temp_perf_write.insert(0, perf_key)
                    perf_write.writerow(temp_perf_write)

    if strict_in_dict:
        with open(strict_file_name, 'w', newline='\n') as strict_file:
            strict_write = csv.writer(strict_file, delimiter=',')
            strict_write.writerow(["Best_Hit_ARO", "ARO", "Cut_Off",  "Sum_Hits", "Model_Type"])
            for strict_key, strict_item in rgi_dict.items():
                if strict_item[1] == "Strict":
                    temp_strict_write = strict_item.copy()
                    temp_strict_write.insert(0, strict_key)
                    strict_write.writerow(temp_strict_write)

    if loose_in_dict:
        with open(loose_file_name, 'w', newline='\n') as loose_file:
            loose_write = csv.writer(loose_file, delimiter=',')
            loose_write.writerow(["Best_Hit_ARO", "ARO", "Cut_Off",  "Sum_Hits", "Model_Type"])
            for loose_key, loose_item in rgi_dict.items():
                if loose_item[1] == "Loose":
                    temp_loose_write = loose_item.copy()
                    temp_loose_write.insert(0, loose_key)
                    loose_write.writerow(temp_loose_write)


if __name__ == '__main__':
    rgi_output(sys.argv[1])
