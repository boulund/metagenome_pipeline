#!/usr/bin/env python
# Find what FASTQ files to combine to produce complete samples representing
# single individuals at single time points
# Fredrik Boulund 2017

from sys import argv, exit
import os
import argparse
from collections import defaultdict
from fnmatch import fnmatch


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("DIR", nargs="+", help="Directory to scan for files.")
    parser.add_argument("--read1", default="*_R1_*", 
            help="Pattern to identify first read of pair ['%(default)s'].")
    parser.add_argument("--read2", default="*_R2_*", 
            help="Pattern to identify second read of pair ['%(default)s'].")
    parser.add_argument("--pattern", default="*.fastq.gz",
            help="Pattern to identify files ['%(default)s'].")
    parser.add_argument("--dest", default="concat_individual_timepoints",
            help="Destination dir ['%(default)s']")
    if len(argv) < 2:
        parser.print_help()
        exit(1)
    return parser.parse_args()


def find_files(directories, pattern):
    """
    Recursively search dirs with a glob pattern yielding paths to matching filenames.

    :param directory:  Path to directory from which to start walking.
    :param pattern:  a glob pattern used select what files to find.
    :return:  yields path to each filename that matches the glob pattern.
    """
    for directory in directories:
        for root, subfolders, files in os.walk(directory, followlinks=True):
            for basename in files:
                if fnmatch(basename, pattern):
                    filename = os.path.join(root, basename)
                    yield root, basename 


def group_files_on_individual_timepoint(filenames):
    grouped_files = defaultdict(list)
    for root, filename in filenames:
        individual_timepoint = filename.split("_", maxsplit=1)[0]
        grouped_files[individual_timepoint].append(os.path.join(root, filename))
    return grouped_files


def print_file_merge_commands(grouped, destination, read_number):
    print("mkdir -p {dest}".format(dest=destination))
    for individual, files in grouped.items():
        fn_string = " ".join(files)
        dest_fn = os.path.join(destination, individual+"_R"+str(read_number)+".fastq.gz")
        if len(files) > 1:
            cmd = "zcat"
            redirect = ">"
            comp = "| gzip "
        else:
            cmd = "cp"
            redirect = ""
            comp = ""
        
        print("{command} {filenames} {comp} {redir} {destination}".format(command=cmd,
                                                             filenames=fn_string,
                                                             comp=comp,
                                                             redir=redirect,
                                                             destination=dest_fn))


if __name__ == "__main__":
    options = parse_args()
    pattern = options.pattern
    read1 = options.read1
    read2 = options.read2
    directory = options.DIR
    destination = options.dest

    files = sorted(list(find_files(directory, pattern)), key=lambda x: x[1])
    ones = [f for f in files if fnmatch(f[1], read1)]
    twos = [f for f in files if fnmatch(f[1], read2)]

    ones_grouped = group_files_on_individual_timepoint(ones)
    twos_grouped = group_files_on_individual_timepoint(twos)

    print_file_merge_commands(ones_grouped, destination, "1")
    print_file_merge_commands(twos_grouped, destination, "2")
