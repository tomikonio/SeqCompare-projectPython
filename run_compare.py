#!/bin/python
import os
import shutil
import subprocess




def set_cmd():
    os.chdir(dir)

    new_dir = "{}/compare{}".format(dir, current_compare)

    if current_compare == 1:
        os.mkdir(new_dir)

    # Todo copies
    shutil.copyfile(src=input_nucl, dst="{}/{}".format(new_dir, input_nucl))
    shutil.copyfile(src=input_protein, dst="{}/{}".format(new_dir, input_protein))

    os.chdir(new_dir)

    makeblastdb_cline = "makeblastdb -in {} -dbtype nucl".format(input_nucl)
    tblastn_cline = "tblastn -query {} -db {} -out compare_output{}.txt".format(input_protein, input_nucl,
                                                                                current_compare)
    next_input = "matching_protein_fasta" if compare_files[input_nucl] == "m" else "not_matching_protein_fasta"

    return makeblastdb_cline, tblastn_cline, next_input



dir = os.getcwd()
# files to be compared with
compare_files = {"fgp2.fasta": "m"}
# compare_files = {"fgp2.fasta": "m", "9802.fasta": "m", "6803.fasta": "nm", "29413.fasta": "nm", "7120.fasta": "nm",
#                  "7942.fasta": "nm"}
# match_type = ("m", "m", "nm", "nm", "nm", "nm")
initial_file = "leptolyngbya.fasta"

total_compares = len(compare_files)
current_compare = 1
input_protein = initial_file

for input_nucl in compare_files:
    makeblastdb_cline, tblastn_cline, next_input = set_cmd()

    subprocess.run(makeblastdb_cline.split())
    subprocess.run(tblastn_cline.split())
    print(input_nucl)
    print(compare_files[input_nucl])
