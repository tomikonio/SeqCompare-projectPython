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


def push_to_match_dict(key):
    match_dict[key] = "1"


def push_to_not_match_dict(key):
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
    FILEHANDLE1.write(columns)
    FILEHANDLE2.write(columns)

    push_fasta_contigs_to_dict()

    for query_result in input_protein:
        query_name = query_result.id
        query_description = query_result.description
        full_name = "{} {}".format(query_name, query_description)

        # output "no hits found" if there is no hits
        if len(query_result.hit_keys) == 0:
            # if this protein already exists as not matched, write it to repeating_protein.txt
            if full_name in not_match_dict:
                FILEHANDLE3.write("{}\t".format(full_name))
                FILEHANDLE3.write("\n")
            else:
                global count_found_no_hits
                count_found_no_hits += 1
                FILEHANDLE2.write("{} {}---===***No Hits Found***===---\n".format(full_name, query_result.seq_len))
                # Todo check if len(query_result) is the same as $result->query_length in perl

                # write to dict of all not matching
                push_to_not_match_dict(full_name)

                # write to not_matching_protein_fasta
                sequence_manipulations(full_name, query_result.seq_len, out_file_not_matching, alphabet)

        # hits found
        else:
            temp_evalue = -1
            temp_hsp = ""
            hit = next(query_result.hits)
            # Todo check if next(query_result.hits) is the same as $result->next_hit in perl
            hsp = next(hit.hsps)
            query_cover_len = hsp.query_span
            # Todo check if hsp.query_span is the same as $hsp->end('query') - $hsp->start('query') in perl

            # Select only the best hsp hit-> # and $highest_score < $hit->bits
            if hsp.evalue <= 1e-040:

                # matched repeating_protein
                if full_name in match_dict:
                    FILEHANDLE3.write("{}\t".format(full_name))
                    FILEHANDLE3.write("\n")
                else:
                    global count_found_hits
                    count_found_hits += 1

                    push_to_match_dict(full_name)
                    sequence_manipulations(full_name, query_result.seq_len, out_file_matching, alphabet)

                    FILEHANDLE1.write(
                        "{} {} {} {} {} {} {} {} {} {}".format(full_name, query_result.seq_len, query_cover_len,
                                                               hit.accession, hit.seq_len, hit.description,
                                                               "hit.segnificance", hsp.bits, hsp.query_frame,
                                                               hsp.query_start, hsp.query_end, hsp.hit_start,
                                                               hsp.hit_end, "$hsp->frac_conserved",
                                                               "$hsp->frac_identical"))
            else:

                if full_name in not_match_dict:
                    FILEHANDLE3.write("{}\t".format(full_name))
                    FILEHANDLE3.write("\n")
                else:
                    global count_found_no_hits
                    count_found_no_hits += 1

                    push_to_not_match_dict(full_name)
                    sequence_manipulations(full_name, query_result.seq_len, out_file_not_matching, alphabet)

                    FILEHANDLE2.write("{}\t{}\n".format(full_name, query_result.seq_len))

                    name_and_length = full_name[7:]  # slice the first 7 characters off the string (why?)
                    # full name in position [7] is the query_result.description, this is a number and thus this number
                    # is made the key in the table dictionary
                    # later, we want output the table dictionary to the FILEHANDELE5 file in an ascending sorted order.
                    parts = name_and_length.split(" ")
                    key = parts[0]
                    value = name_and_length

                    if key not in table:
                        table[key] = []
                    table[key].append(name_and_length)



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

columns = "Query_Name Query_Length Query_Cover_Length Cover_Perc Acc_No Length Desc E_Value Bit_Score Frame QStart QEnd Hit_Start Hit_End Positives Identical"

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

try:
    with open("sorted_not_matching_protein.txt", "w"):
        some_code()
except IOError as e:
    print("Operation failed: {}".format(e.strerror))
