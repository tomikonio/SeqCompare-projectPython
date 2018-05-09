from Bio import SeqIO
from Bio import SearchIO
from flask_table import Table, Col
from beautifultable import BeautifulTable
import sys

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
