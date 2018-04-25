#!/bin/python


# files to be compared with
compare_files = {"fgp2.fasta": "m", "9802.fasta": "m", "6803.fasta": "nm", "29413.fasta": "nm", "7120.fasta": "nm",
                 "7942.fasta": "nm"}
# match_type = ("m", "m", "nm", "nm", "nm", "nm")
initial_file = "leptolyngbya.fasta"

total_compares = len(compare_files)
current_compare = 1
input_protein = initial_file

for input_nucl in compare_files:
    print(input_nucl)