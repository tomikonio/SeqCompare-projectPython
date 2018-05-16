from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from flask_table import Table, Col
from beautifultable import BeautifulTable
import sys



def push_fasta_contigs_to_dict():
    in_file = SeqIO.parse(compare_file,"fasta")

    for result in in_file:
        query_name = "{} ".format(result.id)
        query_description = result.description
        full_name = "{}{}".format(query_name,query_description)

        # Todo full_name regex things and assignment to protein hash


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


protein_dict = {}
match_dict = {}
not_match_dict = {}

table = {}

table = BeautifulTable()
table.column_headers = ["Query_Name", "Query_Length", "Query_Cover_Length", "Cover_Perc Acc_No", "Length", "Desc",
                        "E_Value", "Bit_Score", "Frame", "QStart", "QEnd", "Hit_Start", "Hit_End",
                        "Positives Identical"]
print(table)

input_file = sys.argv[1]
compare_file = sys.argv[2]

# input_handle = open(input_file)
# for query_result in SearchIO.parse(input_handle, format='blast-xml'):
#     # TODO better use the the BLAST libraries in order to parse the blast outputs
#     print()
