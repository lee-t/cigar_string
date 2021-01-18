import sys
import re


# function for reading input file 1. The file is read into a 3D dictionary
# 3D dictionary records the chr start coordinate and cigar string for each transcript-chr alignment
def get_alignment_dict(file_1_name):
    alignment_dict = {}
    file_1 = open(file_1_name, "r")
    for line in file_1:
        line = line.rstrip()
        linearr = re.split(r'\s+', line)

        # check input file 1 format
        check_input_file_1_line_format(linearr)

        trans_id = linearr[0]
        chr_id = linearr[1]
        start_coord = linearr[2]
        cigar = linearr[3]

        alignment_dict[trans_id] = {}
        alignment_dict[trans_id][chr_id] = {}
        alignment_dict[trans_id][chr_id]["start_coord"] = int(start_coord)
        alignment_dict[trans_id][chr_id]["cigar"] = process_cigar_string(cigar)

    file_1.close()
    return alignment_dict


# checks if the current line in the input file 1 has at least 4 fields and the third field contains a integer
# also checks if the cigar string contains valid characters
def check_input_file_1_line_format(inputfile_1_line_arr):
    # checks for number of columns
    if len(inputfile_1_line_arr) < 4:
        print("Error! Inadequate number of fields in input file 1")
        sys.exit()

    # checks if the third column contains an integer
    if not inputfile_1_line_arr[2].isdigit():
        print("Error! No reference coordinate detected")
        sys.exit()

    if not (re.search(r'\d+[MIDNSHP=X]', inputfile_1_line_arr[3], re.I)):
        print("Error! Invalid CIGAR string")
        sys.exit()


# this function takes the cigar string and breaks it into a list of lists
# each element of the outer list is a list containing 2 elements (integer and the cigar character)
def process_cigar_string(cigarstring):
    cigararr = list()
    cigar_num = ''
    for cigar_char in cigarstring:
        if not (cigar_char.isalpha() or cigar_char == '='):
            cigar_num = cigar_num + cigar_char
        else:
            if cigar_num == '':
                print("Error! Wrong CIGAR string format")
                sys.exit()
            cigar_num = int(cigar_num)
            cigararr.append([cigar_num, cigar_char])
            cigar_num = ''
    return cigararr


# function for establishing coordinate correspondence between chr and trans using the 2D cigar arr
# internally calls the "process_cigar_arr" function that returns a dictionary containing coordinate correspondence
def get_coordinate_correspondence(alignment_dict):
    coord_corr = {}
    for tr in alignment_dict:
        coord_corr[tr] = {}
        for ch in alignment_dict[tr]:
            coord_corr[tr][ch] = {}
            chr_start_coord = alignment_dict[tr][ch]["start_coord"]
            cigar_arr = alignment_dict[tr][ch]["cigar"]
            corr_dict = process_cigar_arr(chr_start_coord, cigar_arr)
            coord_corr[tr][ch] = corr_dict
    return coord_corr


# function that returns a dictionary containing coordinate correspondence the key is the trans coord and value is chr

# coord assuming that every trans coord maps uniquely to a chr coord also assumes that cigar contains only 3 chars (
# M,D,I). This function can be easily modified to accommodate other CIGAR characters
def process_cigar_arr(chr_start_coord, cigar_arr):
    corr_dict = {}
    tr_counter = 0
    chr_counter = chr_start_coord
    for cigar_entry in cigar_arr:
        cigar_num = cigar_entry[0]
        cigar_char = cigar_entry[1]

        # alignment operations defined on page 8 of SAM/BAM format specification file
        # these operations consume both query and reference
        if re.match(r'[M=X]', cigar_char, re.I):
            for counter in range(cigar_num):
                corr_dict[tr_counter] = chr_counter
                tr_counter += 1
                chr_counter += 1

        # these operations consume the reference but not the query
        elif re.match(r'[DN]', cigar_char, re.I):
            for counter in range(cigar_num):
                chr_counter += 1

        # these operations consume the query but not the reference
        elif re.match(r'[IS]', cigar_char, re.I):
            for counter in range(cigar_num):
                tr_counter += 1

    return corr_dict


# reads Input file 2 and the 3D coord correspondence dictionary from the "get_coordinate_correspondence" function
# writes the results in the output file which is created in the working directory
def process_input_file_2(file_2_name, coord_corr):
    file_2 = open(file_2_name, "r")
    output_file = open("Output.txt", "w")
    for line in file_2:
        line = line.rstrip()
        linearr = re.split(r'\s+', line)

        # check input file 2 format
        check_input_file_2_line_format(linearr)

        tr_id = linearr[0]
        tr_coord = linearr[1]

        # checks if the transcript id is known
        if not (tr_id in coord_corr):
            print("Warning! Transcript", tr_id, "is not present in Input file 1")
            continue

        # "coord_corr[tr_id][ch_entry][int(tr_coord)])" gets the corresponding chr coord for the given trans coord
        for ch_entry in coord_corr[tr_id]:
            # checks if transcript coordinate is defined
            if not (int(tr_coord) in coord_corr[tr_id][ch_entry]):
                print("Warning! Transcript coordinate", tr_coord, "for transcript", tr_id, "out of alignment range")
                continue

            output_file.write(tr_id + "\t" + tr_coord + "\t" + ch_entry + "\t" + str(
                coord_corr[tr_id][ch_entry][int(tr_coord)]) + "\n")

    file_2.close()
    output_file.close()


# function for checking input file 2 format
def check_input_file_2_line_format(inputfile_2_line_arr):
    # checks for number of columns
    if len(inputfile_2_line_arr) < 2:
        print("Error! Inadequate number of fields in input file 2")
        sys.exit()

    # checks if the second column contains an integer
    if not inputfile_2_line_arr[1].isdigit():
        print("Error! No transcript coordinate detected")
        sys.exit()


# function that executes the entire workflow
def execute_analysis(inputfile_1_name, inputfile_2_name):
    tr_chr_alignment = get_alignment_dict(inputfile_1_name)
    tr_chr_corr = get_coordinate_correspondence(tr_chr_alignment)
    process_input_file_2(inputfile_2_name, tr_chr_corr)
if __name__ == "__main__":
    pass

execute_analysis(sys.argv[1], sys.argv[2])
