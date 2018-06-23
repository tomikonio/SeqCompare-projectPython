#!/bin/python
import os
import shutil
import subprocess
import sys
import json
from collections import OrderedDict
import sort
import compare

"""
This file is the entry point of the script.
It is responsible for managing the folder creation and the execution of the comparisons.
Can be imported as a module, using the start() method from the importing file.
Can be used from the command line:
First argument is the name of the primary file
Second argument is a JSON string of type '{"name": "m", ...}' - where name - secondary filename, m/nm  - match type
Third argument is a full path to the fasta files folder
"""


def set_cmd(current_dir, current_compare, input_protein, input_nucl, matching_type):
    """
    Construct the BLAST commands, generate the next input type (match/not match)

    :param current_dir: the path to the working directory -- the main directory where the files are located
    :param current_compare: int, the number of the compare, starts with 1
    :param input_protein: the primary sequence FASTA file
    :param input_nucl: FASTA file, the secondary nucleotide sequence file
    :param matching_type: string, represents if the input_nucl file`s results are to be match/not match
    :return: makeblastdb_cline, tblastn_cline - two blast commands, next_input - the next input type
    """
    os.chdir(current_dir)

    new_dir = "{}/compare{}".format(current_dir, current_compare)

    # copy the primary and secondary files to the new folder
    if current_compare == 1:
        os.mkdir(new_dir)
        shutil.copyfile(src=input_protein, dst="{}/{}".format(new_dir, input_protein))

    shutil.copyfile(src=input_nucl, dst="{}/{}".format(new_dir, input_nucl))

    os.chdir(new_dir)

    # construct the BLAST commands to be executed
    makeblastdb_cline = "makeblastdb -in {} -dbtype nucl".format(input_nucl)
    tblastn_cline = "tblastn -query {} -db {} -out compare_output{}.xml -outfmt 5".format(input_protein, input_nucl,
                                                                                current_compare)
    next_input = "matching_protein_fasta" if matching_type == "m" else "not_matching_protein_fasta"

    return makeblastdb_cline, tblastn_cline, next_input


def count(next_input, input_nucl, FILEHANDLE1):
    """
    Counts the number of sequences in the output file (the next input file)
    :param next_input: the type of the file to be used as the next primary file
    :param input_nucl: secondary file name
    :param FILEHANDLE1: count.txt file handle
    :return:
    """
    counter = 0
    try:
        with open(next_input, "r") as FILEHANDLE2:
            for line in FILEHANDLE2:
                if line[0] == ">":
                    counter += 1
        FILEHANDLE1.write("{}/{}:{}\n".format(input_nucl, next_input, counter))
    except IOError as e:
        print("Operation failed: {}".format(e.strerror))


def start(initial_file, compare_files, current_dir):
    """
    Start the script, responsible for creating folders and executing the compare file
    :param initial_file: the first primary file
    :param compare_files: an OrderedDict of secondary files and match types
    :param current_dir: a full path to the directory of the fasta files
    :return:
    """

    total_compares = len(compare_files)
    current_compare = 1
    input_protein = initial_file

    with open("{}/counts.txt".format(current_dir), "w") as FILEHANDLE1:
        print("creating counts file...")

        for input_nucl in compare_files:
            makeblastdb_cline, tblastn_cline, next_input = set_cmd(current_dir, current_compare, input_protein, input_nucl, compare_files[input_nucl])

            subprocess.run(makeblastdb_cline.split())
            subprocess.run(tblastn_cline.split())

            compare.start("compare_output{}.xml".format(current_compare), input_protein)
            sort.sort(next_input, "sorted_{}.txt".format(next_input))


            count(next_input, input_nucl, FILEHANDLE1)

            next_compare = current_compare + 1
            if next_compare <= total_compares:
                next_file = "{}_output{}.fasta".format(next_input, current_compare)     # the next primary file

                os.mkdir("{}/compare{}".format(current_dir, next_compare))

                shutil.copyfile(next_input, "{}/compare{}/{}".format(current_dir, next_compare, next_file))
                input_protein = next_file
                current_compare += 1
                print(input_nucl)
                print(compare_files[input_nucl])


if __name__ == '__main__':
    primary_file = sys.argv[1]
    json_dict = sys.argv[2]
    folder_path = sys.argv[3]

    parsed_dict = json.loads(json_dict, object_pairs_hook=OrderedDict)
    start(primary_file, parsed_dict, folder_path)
