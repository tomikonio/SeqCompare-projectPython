from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import csv
import sys

"""
This file is the main script file in terms of functionality,
responsible for making the actual comparison between the primary file
and the BLAST output.
Can be used as a module - activating the start() method from the importing file.
Can be used separetly from the command line.
"""


def push_fasta_contigs_to_dict(compare_file, protein_dict):
    """
    Populate the protein dictionary with sequence entries
    :param compare_fle: name of the current primary file
    :param protein_dict: a dictionary containing key-value pairs of (name, sequence) 
    :return:
    """
    in_file = SeqIO.parse(compare_file, "fasta")

    for record in in_file:
        query_description = record.description
        full_name = "{}".format(query_description)
        # Substitute multiple spaces with one space and remove trailing spaces.
        full_name = ' '.join(full_name.split())

        if full_name not in protein_dict:
            protein_dict[full_name] = record.seq


def push_to_match_dict(key, match_dict):
"""
Push key value pairs of (name, 1) to the match_dict
:param key: key of a dictionary - a full name of a sequence
:param match_dict: dictionary containing enteries that signify that they have been matched with the primary file
:return:
"""
    match_dict[key] = "1"


def push_to_not_match_dict(key, not_match_dict):
"""
Push key value pairs of (name, 0) to the not_match_dict
:param key: key of a dictionary - a full name of a sequence
:param not_match_dict: dictionary containing enteries that signify that they have not been matched with the primary file
:return:
"""
    not_match_dict[key] = "0"


def sequence_manipulations(key, query_length, output_file, alphabet, protein_dict):
    """
    Constructs a SeqRecord object and writes it the given output_file - only if the key is present in the protein_dict
    :param key: key of a dictionary - a full name of a sequence
    :param query_length: length of a query
    :param output_file: the file to write to - matching file/not matching file
    :param alphabet: alphabet of the sequence
    :param protein_dict: dictionary containing full_name as keys and sequences as values
    :return:
    """
    query_length_string = "length={}".format(query_length)

    # construct a SeqRecord and write it to the given file
    if key in protein_dict:
        seqrecord_obj = SeqRecord(Seq(str(protein_dict[key]), alphabet))
        seqrecord_obj.description = key
        seqrecord_obj.id = ""
        SeqIO.write(seqrecord_obj, output_file, "fasta")


def make_compare(compare_file, table, input_protein, out_file_matching, out_file_not_matching, FILEHANDLE1, FILEHANDLE2,
                 FILEHANDLE3):

    protein_dict = {}
    match_dict = {}
    not_match_dict = {}

    table_headers = ["Query_Name", "Query_Length", "Query_Cover_Length", "Acc_No", "Length", "Desc",
                     "E_Value", "Bit_Score", "Frame", "QStart", "QEnd", "Hit_Start", "Hit_End",
                     "Positives","Identical"]

    csv_writer_matching = csv.writer(FILEHANDLE1)
    csv_writer_not_matching = csv.writer(FILEHANDLE2)

    csv_writer_matching.writerow(table_headers)
    csv_writer_not_matching.writerow(table_headers)

    count_found_hits = 0
    count_found_no_hits = 0

    alphabet = "protein"

    push_fasta_contigs_to_dict(compare_file, protein_dict)

    for blast_record in input_protein:

        full_name = blast_record.query

        # output "no hits found" if there is no hits
        if len(blast_record.alignments) == 0:
            # if this protein already exists as not matched, write it to repeating_protein.txt
            if full_name in not_match_dict:
                FILEHANDLE3.write("{}\t".format(full_name))
                FILEHANDLE3.write("\n")
            else:

                count_found_no_hits += 1
                csv_writer_not_matching.writerow([full_name, blast_record.query_letters, "---===***No Hits Found***===---"])

                # write to dict of all not matching
                push_to_not_match_dict(full_name, not_match_dict)

                # write to not_matching_protein_fasta
                sequence_manipulations(full_name, blast_record.query_letters, out_file_not_matching, alphabet,
                                       protein_dict)

        # hits found
        else:

            hit = blast_record.alignments[0]
            hsp = hit.hsps[0]
            query_cover_len = hsp.query_end - hsp.query_start

            # Select only the best hsp hit
            if hsp.expect <= 1e-040:

                # matched repeating_protein - write it to repeating_protein.txt
                if full_name in match_dict:
                    FILEHANDLE3.write("{}\t".format(full_name))
                    FILEHANDLE3.write("\n")
                else:

                    count_found_hits += 1

                    push_to_match_dict(full_name, match_dict)
                    sequence_manipulations(full_name, blast_record.query_letters, out_file_matching, alphabet,
                                           protein_dict)


                    """
                    Comparing to BioPerl:
                    full_name - $result->query_name . " ",$result->query_description
                    blast_record.query_letters - $result->query_length
                    query_cover_len = hsp.query_end - hsp.query_start - $hsp->end('query') - $hsp->start('query')
                    hit.accession - $hit->accession
                    hit.length - $hit->length
                    hit.title - $hit->description   ## hit.title actualy gives more information
                    blast_record.descriptions[0].e - $hit->significance
                    hsp.bits - $hsp->bits
                    hsp.frame - $hsp->query->frame      ## hsp.frame outputs a tuple (num, num), while $hsp->query->frame outputs a number
                    hsp.query_start - $hsp->start('query')
                    hsp.query_end - $hsp->end('query')
                    hsp.sbjct_start - $hsp->start('hit')    ## hsp.sbjct_start - The start residue for the sbjct sequence, is not the same as hit start in perl.
                    hsp.sbjct_end - $hsp->end('hit')        ## same as the above
                    hsp.positives - $hsp->frac_conserved    ## hsp.positives outputs the number of positives, while $hsp->frac_conserved outputs the fraction of positive hits.
                    hsp.identities - $hsp->frac_identical   ## same as the above

                    """
                    # blast_record.descriptions[0].e is the expected value of the hit, supposed to correspond to $hit->significance in perl

                    row = [full_name, blast_record.query_letters,
                                   query_cover_len,
                                   hit.accession, hit.length,
                                   hit.title,
                                   blast_record.descriptions[0].e, hsp.bits,
                                   hsp.frame,
                                   hsp.query_start, hsp.query_end,
                                   hsp.sbjct_start,
                                   hsp.sbjct_end, hsp.positives,
                                   hsp.identities]
                    
                    # write to matching protein.xslx
                    csv_writer_matching.writerow(row)

            else:
                # write to repeating_protein.txt
                if full_name in not_match_dict:
                    FILEHANDLE3.write("{}\t".format(full_name))
                    FILEHANDLE3.write("\n")
                else:
                    # consider the hit as a miss beacuse it is above 1e-040
                    count_found_no_hits += 1

                    push_to_not_match_dict(full_name, not_match_dict)
                    sequence_manipulations(full_name, blast_record.query_letters, out_file_not_matching, alphabet,
                                           protein_dict)

                    csv_writer_not_matching.writerow([full_name, blast_record.query_letters])

                    name_and_length = "{}\t {}\n".format(full_name, blast_record.query_letters)
                    fragment = name_and_length[7:]  # slice the first 7 characters off the string

                    # full name in position [7] is the start of the blast_record.query_letters, this is a number and thus this number
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

    table = {}

    blast_type = "blast-xml"  # blast-xml
    # input_protein = SearchIO.parse(input_file, format=blast_type)

    with open(input_file, "r") as input_protein_handle, open("matching_protein_fasta", "w") as out_file_matching, open(
            "not_matching_protein_fasta",
            "w") as out_file_not_matching:
        print("making files....")
        # out_file_matching = "matching_protein_fasta"
        # out_file_not_matching = "not_matching_protein_fasta"

        input_protein = NCBIXML.parse(input_protein_handle)

        try:
            with open("matching_protein.xlsx", "w") as FILEHANDLE1, open("not_matching_protein.xlsx",
                                                                        "w") as FILEHANDLE2, open(
                "repeating_protein.txt", "w") as FILEHANDLE3:

                make_compare(compare_file, table, input_protein, out_file_matching, out_file_not_matching, FILEHANDLE1,
                             FILEHANDLE2, FILEHANDLE3)

        except IOError as e:
            print("Operation failed: {}".format(e.strerror))

        try:
            with open("sorted_not_matching_protein.txt", "w") as FILEHANDLE5:
                make_sorted_not_matching_file(table, FILEHANDLE5)
        except IOError as e:
            print("Operation failed: {}".format(e.strerror))




if __name__ == "__main__":
    input_file = sys.argv[1]
    compare_file = sys.argv[2]
    start(input_file, compare_file)
