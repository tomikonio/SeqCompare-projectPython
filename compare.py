from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from flask_table import Table, Col
# from beautifultable import BeautifulTable
import sys


def push_fasta_contigs_to_dict(compare_file, protein_dict):
    """
    Populate the protein dictionary with sequence entries
    :return:
    """
    in_file = SeqIO.parse(compare_file, "fasta")

    for record in in_file:
        # query_name = "{} ".format(record.id)
        # print("record decription is: "+ record.description)
        query_description = record.description
        # full_name = "{} {}".format(query_name, query_description)
        full_name = "{}".format(query_description)
        # print("Full name is " + full_name)
        # Substitute multiple spaces with one space and remove trailing spaces.
        full_name = ' '.join(full_name.split())

        if full_name not in protein_dict:
            protein_dict[full_name] = record.seq
    # print(" ")


def push_to_match_dict(key, match_dict):
    match_dict[key] = "1"


def push_to_not_match_dict(key, not_match_dict):
    not_match_dict[key] = "0"


def sequence_manipulations(key, query_length, output_file, alphabet, protein_dict):
    """

    :param key: key of a dictionary
    :param query_length: length of a query
    :param output_file: the file to write to
    :param alphabet: alphabet of the sequence
    :param protein_dict: dictionary containing full_name as keys and sequences as values
    :return:
    """
    query_length_string = "length={}".format(query_length)
    # dict_val = protein_dict[key]

    # print("The key is is: "+ key)
    if key in protein_dict:
        seqrecord_obj = SeqRecord(Seq(str(protein_dict[key]), alphabet))
        seqrecord_obj.description = key
        seqrecord_obj.id = ""
        SeqIO.write(seqrecord_obj, output_file, "fasta")


def make_compare(compare_file, protein_dict, match_dict, not_match_dict, table, columns,
                 input_protein, count_found_hits, count_found_no_hits, alphabet, out_file_matching,
                 out_file_not_matching, FILEHANDLE1, FILEHANDLE2, FILEHANDLE3):
    # global count_found_no_hits
    # global count_found_hits

    FILEHANDLE1.write(columns)
    FILEHANDLE2.write(columns)

    push_fasta_contigs_to_dict(compare_file, protein_dict)

    for query_result in input_protein:
        query_name = query_result.id
        query_description = query_result.description
        full_name = "{} {}".format(query_name, query_description)

        # output "no hits found" if there is no hits
        if len(query_result.hit_keys) == 0:
            # if this protein already exists as not matched, write it to repeating_protein.txt
            if full_name in not_match_dict:
                print("repeating......")
                FILEHANDLE3.write("{}\t".format(full_name))
                FILEHANDLE3.write("\n")
            else:

                count_found_no_hits += 1
                FILEHANDLE2.write("{} {}---===***No Hits Found***===---\n".format(full_name, query_result.seq_len))
                # Todo check if len(query_result) is the same as $result->query_length in perl

                # write to dict of all not matching
                push_to_not_match_dict(full_name, not_match_dict)

                # write to not_matching_protein_fasta
                sequence_manipulations(full_name, query_result.seq_len, out_file_not_matching, alphabet, protein_dict)

        # hits found
        else:
            temp_evalue = -1
            temp_hsp = ""
            hit = query_result.hits[0]
            # Todo check if next(query_result.hits) is the same as $result->next_hit in perl
            hsp = hit.hsps[0]
            query_cover_len = hsp.query_span
            # Todo check if hsp.query_span is the same as $hsp->end('query') - $hsp->start('query') in perl

            # Select only the best hsp hit-> # and $highest_score < $hit->bits
            if hsp.evalue <= 1e-040:

                # matched repeating_protein
                if full_name in match_dict:
                    print("repeating......")
                    FILEHANDLE3.write("{}\t".format(full_name))
                    FILEHANDLE3.write("\n")
                else:

                    count_found_hits += 1

                    push_to_match_dict(full_name, match_dict)
                    sequence_manipulations(full_name, query_result.seq_len, out_file_matching, alphabet, protein_dict)

                    FILEHANDLE1.write(
                        "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(full_name, query_result.seq_len,
                                                                                query_cover_len,
                                                                                hit.accession, hit.seq_len,
                                                                                hit.description,
                                                                                "hit.segnificance", hsp.bitscore,
                                                                                hsp.query_frame,
                                                                                hsp.query_start, hsp.query_end,
                                                                                hsp.hit_start,
                                                                                hsp.hit_end, hsp.pos_num,
                                                                                hsp.ident_num))
            else:

                if full_name in not_match_dict:
                    print("repeating......")
                    FILEHANDLE3.write("{}\t".format(full_name))
                    FILEHANDLE3.write("\n")
                else:

                    count_found_no_hits += 1

                    push_to_not_match_dict(full_name, not_match_dict)
                    sequence_manipulations(full_name, query_result.seq_len, out_file_not_matching, alphabet, protein_dict)

                    FILEHANDLE2.write("{}\t{}\n".format(full_name, query_result.seq_len))

                    name_and_length = "{}\t {}\n".format(full_name, query_result.seq_len)
                    fragment = name_and_length[7:]  # slice the first 7 characters off the string (why?)

                    # full name in position [7] is the start of the query_result.description, this is a number and thus this number
                    # is made the key in the table dictionary
                    # later, we want output the table dictionary to the FILEHANDELE5 file in an ascending sorted order.

                    parts = fragment.split(" ")
                    key = parts[0]
                    value = name_and_length

                    if key not in table:
                        table[key] = []
                    table[key].append(value)


def make_sorted_not_matching_file(table, FILEHANDLE5):
    """
    Outputs the values held in the table dict to the make_sorted_not_matching_protein.txt file
    Each line is represents a sequence
    Because one dict entry may hold several sequence names, each dict entry is iterated over
    :return:
    """
    key_list = list(table.keys())
    key_list.sort(key=int)  # sort the dict keys by ascending order

    for key in key_list:
        for sequence in table[key]:
            FILEHANDLE5.write("{}\n".format(sequence))


def start(input_file, compare_file):
    print("Started the compare file")
    protein_dict = {}
    match_dict = {}
    not_match_dict = {}

    table = {}

    columns = "Query_Name Query_Length Query_Cover_Length Cover_Perc Acc_No Length Desc E_Value Bit_Score Frame QStart QEnd Hit_Start Hit_End Positives Identical\n"
    blast_type = "blast-xml"  # blast-xml
    input_protein = SearchIO.parse(input_file, format=blast_type)
    progress = 0
    count_found_hits = 0
    count_found_no_hits = 0

    alphabet = "protein"

    with open("matching_protein_fasta", "w") as out_file_matching, open("not_matching_protein_fasta",
                                                                        "w") as out_file_not_matching:
        print("making files....")
        # out_file_matching = "matching_protein_fasta"
        # out_file_not_matching = "not_matching_protein_fasta"

        try:
            with open("matching_protein.txt", "w") as FILEHANDLE1, open("not_matching_protein.txt",
                            "w") as FILEHANDLE2, open("repeating_protein.txt", "w") as FILEHANDLE3:

                make_compare(compare_file, protein_dict, match_dict, not_match_dict, table, columns,
                             input_protein, count_found_hits, count_found_no_hits, alphabet, out_file_matching,
                             out_file_not_matching, FILEHANDLE1, FILEHANDLE2, FILEHANDLE3)

        except IOError as e:
            print("Operation failed: {}".format(e.strerror))

        try:
            with open("sorted_not_matching_protein.txt", "w") as FILEHANDLE5:
                make_sorted_not_matching_file(table, FILEHANDLE5)
        except IOError as e:
            print("Operation failed: {}".format(e.strerror))


# protein_dict = {}
# match_dict = {}
# not_match_dict = {}
#
# table = {}
# print("Started the compare file")
# # table = BeautifulTable()
# # table.column_headers = ["Query_Name", "Query_Length", "Query_Cover_Length", "Cover_Perc Acc_No", "Length", "Desc",
# #                         "E_Value", "Bit_Score", "Frame", "QStart", "QEnd", "Hit_Start", "Hit_End",
# #                         "Positives Identical"]
# # print(table)
# #
# input_file = sys.argv[1]
# compare_file = sys.argv[2]
#
# # input_handle = open(input_file)
# # for query_result in SearchIO.parse(input_handle, format='blast-xml'):
# #     # TODO better use the the BLAST libraries in order to parse the blast outputs
# #     print()
#
# columns = "Query_Name Query_Length Query_Cover_Length Cover_Perc Acc_No Length Desc E_Value Bit_Score Frame QStart QEnd Hit_Start Hit_End Positives Identical\n"
# format = "blast-xml"    #blast-xml
# input_protein = SearchIO.parse(input_file, format=format)
# progress = 0
# count_found_hits = 0
# count_found_no_hits = 0
#
# alphabet = "protein"
#
# with open("matching_protein_fasta", "w") as out_file_matching, open("not_matching_protein_fasta", "w") as out_file_not_matching:
#     print("making files....")
# # out_file_matching = "matching_protein_fasta"
# # out_file_not_matching = "not_matching_protein_fasta"
#
#     try:
#         with open("matching_protein.txt", "w") as FILEHANDLE1, open("not_matching_protein.txt", "w") as FILEHANDLE2, open("repeating_protein.txt", "w") as FILEHANDLE3:
#
#             make_compare()
#
#     except IOError as e:
#         print("Operation failed: {}".format(e.strerror))
#
#     try:
#         with open("sorted_not_matching_protein.txt", "w") as FILEHANDLE5:
#             make_sorted_not_matching_file()
#     except IOError as e:
#         print("Operation failed: {}".format(e.strerror))


if __name__ == "__main__":
    input_file = sys.argv[1]
    compare_file = sys.argv[2]
    start(input_file, compare_file)
