import sys


def sort(input_file, output_file):
    table = {}
    with open(input_file, "r") as FILEHANDLE1:

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

    with open(output_file, "w") as FILEHANDLE2:
        key_list = list(table.keys())
        key_list.sort(key=int)  # sort the keys in ascending order
        for key in key_list:
            for sequence in table[key]:
                FILEHANDLE2.write("{}\n".format(sequence))


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sort(input_file, output_file)
