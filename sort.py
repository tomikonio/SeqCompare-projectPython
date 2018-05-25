import sys

input_file = sys.argv[1]

with open(input_file, "r") as FILEHANDLE1:
    table = {}

    for line in FILEHANDLE1.readlines():
        chars = line.split(" ")

        if chars[0] == ">":
            fragment = line[7:]
            parts = fragment.split(" ")
            key = parts[0]
            value = line[1:]

            if key not in table:
                table[key] = []

            table[key].append(value)

    key_list = list(table.keys())
    key_list.sort(key=int)  # sort the keys in ascending order
    for key in key_list:
        for sequence in table[key]:
            FILEHANDLE1.write("{}\n".format(sequence))
