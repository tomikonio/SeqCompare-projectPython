from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from flask_table import Table, Col
from beautifultable import BeautifulTable
import sys


def push_fasta_contigs_to_dict():
    """
    Populate the protein dictionary with sequence entries
    :return:
    """
    in_file = SeqIO.parse(compare_file, "fasta")

    for record in in_file:
        query_name = "{} ".format(record.id)
        query_description = record.description
        full_name = "{}{}".format(query_name, query_description)

        # Substitute multiple spaces with one space and remove trailing spaces.
        full_name = ' '.join(full_name.split())

        if full_name not in protein_dict:
            protein_dict[full_name] = record.seq


def push_to_match_hash(key):
    match_dict[key] = "1"


def push_to_not_match_hasg(key):
    not_match_dict[key] = "0"


def sequence_manipulations(key, query_length, output_file, alphabet):
    """

    :param key: key of a dictionary
    :param query_length: length of a query
    :param output_file: the file to write to
    :param alphabet: alphabet of the sequence
    :return:
    """
    query_length_string = "length={}".format(query_length)
    dict_val = protein_dict[key]

    if protein_dict[key]:
        seq_obj = Seq(dict_val, alphabet)
        SeqIO.write(seq_obj, output_file, "fasta")



def make_compare():
    print()


protein_dict = {}
match_dict = {}
not_match_dict = {}

table = {}

# table = BeautifulTable()
# table.column_headers = ["Query_Name", "Query_Length", "Query_Cover_Length", "Cover_Perc Acc_No", "Length", "Desc",
#                         "E_Value", "Bit_Score", "Frame", "QStart", "QEnd", "Hit_Start", "Hit_End",
#                         "Positives Identical"]
# print(table)
#
input_file = sys.argv[1]
compare_file = sys.argv[2]

# input_handle = open(input_file)
# for query_result in SearchIO.parse(input_handle, format='blast-xml'):
#     # TODO better use the the BLAST libraries in order to parse the blast outputs
#     print()

input_protein = SearchIO.parse(input_file, format="blast-xml")
progress = 0
count_found_hits = 0
count_found_no_hits = 0

alphabet = "protein"
out_file_matching = SeqIO.parse("matching_protein_fasta", "fasta")
out_file_not_matching = SeqIO.parse("not_matching_protein_fasta", "fasta")

try:
    with open("matching_protein.txt", "w") as FILEHANDLE1, open("not_matching_protein.txt", "w") as FILEHANDLE2, open(
            "sorted_not_matching_protein.txt", "w") as FILEHANDLE3, open("repeating_protein.txt", "w") as FILEHANDLE4:

        make_compare()

except IOError as e:
    print("Operation failed: {}".format(e.strerror))
