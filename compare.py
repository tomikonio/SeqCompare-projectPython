from Bio import SeqIO
from Bio import SearchIO
import sys

input_file = sys.argv[1]
compare_file = sys.argv[2]

input_handle = open(input_file)
for query_result in SearchIO.parse(input_handle, format='blast-xml'):
    # TODO how to
    print()
