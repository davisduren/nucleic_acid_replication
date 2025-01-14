import matplotlib.pyplot as plt
import seqfold
import os
from functions import *
import time



#calculate the average length of the sequences returned by the main program
#will also find the longest sequence generated throughout any iteration and return it 
def average_sequence_length_and_longest_seq(filename):
    opened = open(filename, "r")
    sequences = opened.read().splitlines()
    sets = []
    current_set = []
    for seq in sequences:
        if seq.startswith("*"): #ignore our notes written to the output file
            continue
        if seq:
            current_set.append(seq)
        else:
            if current_set:
                sets.append(current_set)
                current_set = []
    if current_set:
        sets.append(current_set)

    avg_lengths = []
    longest_sequence = ""
    for s in sets:
        total_length = 0
        for seq in s:
            total_length += len(seq)
            if len(seq) > len(longest_sequence):
                longest_sequence = seq
        avg_length = total_length / len(s)
        avg_lengths.append(avg_length)

    return avg_lengths, "LONGEST SEQUENCE FOUND: " + longest_sequence


#print(average_sequence_length_and_longest_seq("test_file_1.txt"))

#histogram attempt1
def generate_sequence_length_histogram(filename):
    opened = open(filename, "r")
    sequences = opened.read().splitlines()

    sequence_lengths = []
    for seq in sequences:
        if seq.startswith("*"):
            continue
        else:
            sequence_lengths.append(len(seq))
            

    plt.hist(sequence_lengths, bins=100)
    plt.yscale("log")
    plt.xlabel('Sequence Length ( log scale )')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Histogram')
    plt.savefig('sequence_length_histogram.jpg')  # Save the histogram as a .jpg file
    plt.show()

    opened.close()
# Usage example
#generate_sequence_length_histogram('test_file_1.txt')





"""attempting to create a function to identify the structured
regions after the simulation has been run - it already uses the seqfold 
during the simulation to reduce probability of breaking in structured
regions, but after the iterations are over, we want to be able to 
see the structured regions that have resulted"""

def find_structured_regions_in_file(file_to_write, file_to_read):
    writing = open(file_to_write, "w")
    """writing is a file used to save the sequences"""
    reading_from = open(file_to_read, "r")
    sequences_and_info = reading_from.read().splitlines()

    sequences = []
    for seq in sequences_and_info:
        if seq.startswith("*") or seq == "":
            continue
        else:
            if len(seq) >6: #length cutoff is an arbitrary value here
                sequences.append(seq)

    strands_with_structure = []
    for long_strand in sequences:
        print(long_strand)
        structured_info = seqfold.fold(long_strand)
        #print(structured_info)

        #beginning the replicant structured_bonds function

        
        for struct in structured_info:
            i, j = struct.ij[0][0], struct.ij[0][1] 
            """this gives us the beginning / ending index of where structured region is - 
            ex. 1 10 means that there is structure between the first and 10th nucleotide
            Do we want to save this information? 

            LIMITATION: This assumes only one structured region per strand - i have not 
            yet found a way to remedy this
            """
            if ("HAIRPIN" or "STACK" or "BULGE") in struct.desc:
                strands_with_structure.append(long_strand)
                #print(i,j)
            

    
    
    
    print(strands_with_structure)
            
        
        
        #writing.write(structures)


#find_structured_regions_in_file('new_test_txt.txt', 'test_file_1.txt')


""" creating a function to create a new folder each time the code is run 
this folder will contain the .txt file(s) and the histogram of each run
this is done so each run can be easily separated and stored so we can determine consitency
between runs"""

"""
def create_and_save_output(directory, folder_name, text_files, histogram,):

    output_folder = os.path.join(directory, folder_name)
    os.mkdir(output_folder)
    with open(os.path.join(output_folder, "sequences.txt"), "w") as file:
        file.write(text_files)

    with open(os.path.join(output_folder, "histogram.jpg"), "wb") as file:
        file.write(histogram)

    pass
"""



current_time_struct = time.localtime()
formatted_time = str(time.strftime("%Y-%m-%d %H:%M", current_time_struct))
def create_blank_text_file(base_directory, folder_name,):
    """
    Creates a new text file with a unique name based on time function was run
    Saves it to the location on the harddrive specified by base_directory
    """

    try:
        folder_name = f"{folder_name}"
        output_folder = os.path.join(base_directory, folder_name)
        os.mkdir(output_folder)
        file_name = f"run_{formatted_time}.txt"
        file_path = os.path.join(output_folder, file_name)
        with open(file_path, "w"):
            pass
  
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return False  

base_directory = "/Users/dduren3/Desktop/nucleic acid replication/nucleic_acid_replication-main/attempt1"
custom_folder_name = formatted_time
#create_text_file(base_directory, custom_folder_name)



