import sys

"""
Sort the designated next primary file in ascending order.
The sorted file will not be used as the next primary file.
This script file can be imported as a module, using the sort() method in the importing file.
Can be used form the command line.
"""


def sort(input_file, output_file):
    """
    Sorts the file to be used as next primary file - the sorted file is not used as the next primary
    :param input_file: next input file name
    :param output_file: name of sorted file
    :return:
    """
    table = {}
    with open(input_file, "r") as FILEHANDLE1:

        for line in FILEHANDLE1.readlines():
            # parse the input and write to the table
            chars = line.split(" ")

            if chars[0] == ">":

                fragment = line[9:]     # was 7 originally, but there is additional space
                parts = fragment.split(" ")
                key = parts[0]
                value = line[1:]

                if key not in table:
                    table[key] = []

                table[key].append(value)    # one key may have multiple entries

    with open(output_file, "w") as FILEHANDLE2:
        key_list = list(table.keys())
        key_list.sort(key=int)  # sort the keys in ascending order
        for key in key_list:
            for sequence in table[key]:
                FILEHANDLE2.write("{}\n".format(sequence.lstrip()))     # lstrip to remove the leading whitespace


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sort(input_file, output_file)
