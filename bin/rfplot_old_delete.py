import matplotlib.pyplot as plt
import csv
import os
import argparse

# Author: Jacob Singer (https://github.com/JacobSinger42)
# retrieves the command-line arguments and flags

parser = argparse.ArgumentParser()

parser.add_argument('--dir', type=str, help='path (relative or absolute) to the directory with the sample files')
parser.add_argument('--nd', action='store_true', help='no display: whether to display the generated graphs, typically used with -s')
parser.add_argument('--s', action='store_true', help='save: whether to save the generated graphs')
parser.add_argument('--sd', type=str, help='save directory: path (relative or absolute) at which to save the graphs if -s is used')

args = parser.parse_args()

# empty initializations that are later sorted with values
x_gene, y_gene = [], []
x_group, y_group = [], []
x_mech, y_mech = [], []
x_class, y_class = [], []
x_type, y_type = [], []

names = []

# gets data from sample files in directory

dir = args.dir if args.dir[-1] == '/' else args.dir + '/'

for fn in os.listdir(dir):

    f = fn.split('.') # subdivides string at . in filename

    if (f[-1] != 'tsv'): # ensures the file is a .tsv
        continue

    if f[0] not in names:
        names.append(f[0]) # adds the sample name to the list

    # checks the feature type and aliases to the corresponding lists
    if (f[-2] == 'gene'):
        x, y = x_gene, y_gene
    elif (f[-2] == 'group'):
        x, y = x_group, y_group
    elif (f[-2] == 'mech'):
        x, y = x_mech, y_mech
    elif (f[-2] == 'class'):
        x, y = x_class, y_class
    elif (f[-2] == 'type'): 
        x, y = x_type, y_type

    # adds a new empty list to the feature type for the sample
    x.append([])
    y.append([])

    # writes the data from the file to the list
    with open(dir + fn) as file:
        tsv = csv.reader(file, delimiter='\t')
        for i in tsv:
            x[-1].append(int(i[0]))
            y[-1].append(int(i[1]))

# displays retrieved data using PyPlot

sps = []
for i in range (5):
    sps.append(plt.subplots()) # creates a list of tuples with the figure and axes objects

# plots the data from each sample list
for i in range (len(names)): 
    sps[0][1].plot(x_gene[i], y_gene[i], label=names[i])
    sps[1][1].plot(x_group[i], y_group[i], label=names[i])
    sps[2][1].plot(x_mech[i], y_mech[i], label=names[i])
    sps[3][1].plot(x_class[i], y_class[i], label=names[i])
    sps[4][1].plot(x_type[i], y_type[i], label=names[i])

# adds additional formatting to the graphs
for i in range(5):
    sps[i][1].set_xlabel('% of data subsampled')
    sps[i][1].set_ylabel('unique features identified')
    sps[i][1].legend(bbox_to_anchor=(1.05,1.0), loc='upper left')
    sps[i][1].set_xlim(left=0)
    sps[i][1].set_ylim(bottom=0)

    sps[i][0].set_facecolor('white')
    sps[i][0].set_tight_layout(True)

# sets titles for the graphs
sps[0][1].set_title('Gene Subsampling Features')
sps[1][1].set_title('Group Subsampling Features')
sps[2][1].set_title('Mech Subsampling Features')
sps[3][1].set_title('Class Subsampling Features')
sps[4][1].set_title('Type Subsampling Features')

# displays graphing windows
if (not args.nd):
    plt.show()

# adds the save directory, if applicable
if (args.sd):
    sd = args.sd if args.sd[-1] == '/' else args.sd + '/'
else:
    sd = ''

# saves the plots
if (args.s):
    sps[0][0].savefig(sd + 'Gene.png')
    sps[1][0].savefig(sd + 'Group.png')
    sps[2][0].savefig(sd + 'Mech.png')
    sps[3][0].savefig(sd + 'Class.png')
    sps[4][0].savefig(sd + 'Type.png')
