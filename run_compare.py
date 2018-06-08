#!/bin/python
import os
import shutil
import subprocess
from Bio.Blast.Applications import NcbitblastnCommandline
import sort
import compare




def set_cmd(current_dir, current_compare, input_protein, input_nucl, compare_files):
    """
    Construct the BLAST commands, generate the next input type (match/not match)
    :return: makeblastdb_cline, tblastn_cline - two blast commands, next_input - the next input type
    """
    os.chdir(current_dir)

    new_dir = "{}/compare{}".format(current_dir, current_compare)

    if current_compare == 1:
        os.mkdir(new_dir)
        shutil.copyfile(src=input_protein, dst="{}/{}".format(new_dir, input_protein))

    # Todo copie compare.pl - maybe not needed here
    # shutil.copyfile(src="compare.py", dst="{}/compare.py".format(new_dir))
    shutil.copyfile(src=input_nucl, dst="{}/{}".format(new_dir, input_nucl))
    #shutil.copyfile(src=input_protein, dst="{}/{}".format(new_dir, input_protein))

    os.chdir(new_dir)

    makeblastdb_cline = "makeblastdb -in {} -dbtype nucl".format(input_nucl)
    tblastn_cline = "tblastn -query {} -db {} -out compare_output{}.xml -outfmt 5".format(input_protein, input_nucl,
                                                                                current_compare)
    next_input = "matching_protein_fasta" if compare_files[input_nucl] == "m" else "not_matching_protein_fasta"

    return makeblastdb_cline, tblastn_cline, next_input

def count(next_input, input_nucl, FILEHANDLE1):
    """
    Counts the number of sequences in the output file (the next input file)
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
    print(current_dir)
    # current_dir = os.getcwd()
    # files to be compared with
    # compare_files = {"fgp2.fasta": "m", "9802.fasta": "m"}
    # compare_files = {"fgp2.fasta": "m", "9802.fasta": "m", "6803.fasta": "nm", "29413.fasta": "nm", "7120.fasta": "nm",
    #                  "7942.fasta": "nm"}
    # match_type = ("m", "m", "nm", "nm", "nm", "nm")
    # initial_file = "leptolyngbya.fasta"

    total_compares = len(compare_files)
    current_compare = 1
    input_protein = initial_file

    with open("counts.txt", "w") as FILEHANDLE1:
        print("creating counts file...")

        for input_nucl in compare_files:
            makeblastdb_cline, tblastn_cline, next_input = set_cmd(current_dir, current_compare, input_protein, input_nucl, compare_files)

            subprocess.run(makeblastdb_cline.split())
            subprocess.run(tblastn_cline.split())

            # Todo here call the compare.pl equivalent
            # subprocess.run(["python", "compare.py", "compare_output{}.xml".format(current_compare), input_protein])
            compare.start("compare_output{}.xml".format(current_compare), input_protein)
            sort.sort(next_input, "sorted_{}.txt".format(next_input))

            # Todo call the count function - what for?

            count(next_input, input_nucl, FILEHANDLE1)

            # if( $currentCompare <= $totalCompares )
            next_compare = current_compare + 1
            if next_compare <= total_compares:
                next_file = "{}_output{}.fasta".format(next_input, current_compare)

                os.mkdir("{}/compare{}".format(dir, next_compare))

                shutil.copyfile(next_input, "{}/compare{}/{}".format(dir, next_compare, next_file))
                input_protein = next_file
                current_compare += 1
                print(input_nucl)
                print(compare_files[input_nucl])





# dir = os.getcwd()
# # files to be compared with
# compare_files = {"fgp2.fasta": "m", "9802.fasta": "m"}
# # compare_files = {"fgp2.fasta": "m", "9802.fasta": "m", "6803.fasta": "nm", "29413.fasta": "nm", "7120.fasta": "nm",
# #                  "7942.fasta": "nm"}
# # match_type = ("m", "m", "nm", "nm", "nm", "nm")
# initial_file = "leptolyngbya.fasta"
#
# total_compares = len(compare_files)
# current_compare = 1
# input_protein = initial_file
#
# with open("counts.txt", "w") as FILEHANDLE1:
#     print("creating counts file...")
#
#     for input_nucl in compare_files:
#         makeblastdb_cline, tblastn_cline, next_input = set_cmd()
#
#         subprocess.run(makeblastdb_cline.split())
#         subprocess.run(tblastn_cline.split())
#
#         # Todo here call the compare.pl equivalent
#         # subprocess.run(["python", "compare.py", "compare_output{}.xml".format(current_compare), input_protein])
#         compare.start("compare_output{}.xml".format(current_compare), input_protein)
#         sort.sort(next_input,"sorted_{}.txt".format(next_input))
#
#         # Todo call the count function - what for?
#
#         count()
#
#         # if( $currentCompare <= $totalCompares )
#         next_compare = current_compare + 1
#         if next_compare <= total_compares:
#             next_file = "{}_output{}.fasta".format(next_input, current_compare)
#
#             os.mkdir("{}/compare{}".format(dir, next_compare))
#
#             shutil.copyfile(next_input, "{}/compare{}/{}".format(dir, next_compare, next_file))
#             input_protein = next_file
#             current_compare += 1
#             print(input_nucl)
#             print(compare_files[input_nucl])
