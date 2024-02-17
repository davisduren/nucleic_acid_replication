import matplotlib.pyplot as plt
import seqfold
import os
from functions import *
import time









"""NOVEL 2/16/24
Have histogram get made each progress_report freq of runs 
to be able to visually see when we hit steady state"""

#MODIFIED histogram function to achieve above goal
def generate_sequence_length_histogram(file_path, folder_name, progress_report_num = None):

    opened = open(file_path, "r") # needs to be error contingency file
    sequences = opened.read().splitlines()

    first_histogram = True
    sequence_lengths = []
    for seq in sequences:
        if seq.startswith("!!"):
            first_histogram = False

    if first_histogram:
        for seq in sequences:
            if seq.startswith("*") or not seq.strip():
                continue
            else:
                sequence_lengths.append(len(seq))

    else:
    
        for seq in sequences:
            if not seq.startswith("!!" + str(progress_report_num)):
                continue
            else:
                sequence_lengths.append(len(seq))
            

    plt.hist(sequence_lengths, bins=100)
    plt.yscale("log")
    plt.xlabel('Sequence Length (nt)')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Histogram' + str(formatted_time))

    try:
        #"./output/"
        output_folder = "./output/" + str(folder_name)
        output_file_path = os.path.join(output_folder, f"{progress_report_num}_sequence_length_histogram_{formatted_time}.jpg")
        plt.savefig(output_file_path)

    except Exception as e:
        print(e)
        print("The saving of histogram has failed - check output folder?")

    opened.close()

generate_sequence_length_histogram("/Users/dduren3/Desktop/nucleic_acid_replication/output/2024-02-16 19:08/Error_Contingency_run_2024-02-16 19:08.txt", "./output/2024-02-16 19/08", 750 )


print("LINE 65 Running")








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


"""NOVEL ERROR RECOVERY SYSTEM
A note about all this:
I wrote this in about 10 hours after an 11 day run failed and lost it all
What i have here works, but it is not very good
First, you must modify the file name of error_contingency file as follows:
ex. its name is Error_Contingency_run_2024-02-09 22/35.txt
You must change it to Error_Contingency_run_2024-02-09 22-35.txt
See the dash at the end? You must do this manually

"""



inverse_mapping = {"A": "1", "G": "2", "C" : "3", "U" : "4"}


def read_error_save_txt_file(file_path):
    """Read a text file and return its contents as an array of strings.

    Parameters
    ----------
    file_path : str
        Path to the text file.

    Returns
    -------
    np.ndarray
        Array of strings containing the lines of the text file.
    """
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if not (line.startswith('*') or line.startswith("!!"))]
    return np.array(lines)


def convert_str_to_int_seq(seq_array_s, reverse_mapping):
    """Convert string nucleic acid sequences to integers.

    Parameters
    ----------
    seq_array_s : np.ndarray
        Array of string nucleic acid strands
    reverse_mapping : dict
        Map nucleobase names (A, G, C, U) with integer numbers (1, 2, 3, 4)
    
    Returns
    -------
    np.ndarray
        Array of integer nucleic acid strands
    """
    # Convert array elements to strings
    seq_array = seq_array_s.astype(str)
    
    # Replace each base name with the corresponding number
    for k, v in reverse_mapping.items():
        seq_array = [s.replace(k, str(v)) for s in seq_array]
    
    return np.array(seq_array, dtype=int)  # Convert array elements back to integers


# Define the mapping from nucleobase names to str numbers
reverse_mapping = {"A": "1", "G": "2", "C": "3", "U": "4"}
file_path = "/Users/dduren3/Desktop/nucleic_acid_replication/output/2024-02-09 22:35/Error_Contingency_run_2024-02-09 22-35.txt"  # Replace 'your_file_path.txt' with the actual path to your file

#place file path here
sequence_letters = read_error_save_txt_file(file_path)
sequence_numbers = convert_str_to_int_seq(sequence_letters, inverse_mapping)
print(sequence_numbers)




def delete_until_special_chars(file_path, save_num_to_find):
    line_we_need = "!!" + str(save_num_to_find)
    with open(file_path, "r") as file:
        lines = file.readlines()


    index = -1
    settings_str = ""
    for i, line in enumerate(lines):
        if line.startswith("*"):
            settings_str += line
        settings_str += "\n"

        if line.startswith(line_we_need):
            index = i
            break
    
    if index != -1:
        with open(file_path, 'w') as file:
                file.writelines(lines[index:])
    else:
        print("something is amiss with delete function")


def error_continuation(it_failed_at, end_iteration_num, progress_report_freq, file_with_backups, sequence_numbers = None):
    #it failed at == when the code failed (iteration 1308 for ex)
    #end_iteration_num == ultimate goal of how many times to run code (ex. 1500)
    #progress_report_freq == number of it last time code was saved before it failed
    #file_with_backups == name of txt file with backups
    #sequence_numbers == variable thats a numpy array of file containing backups, converted to nums


    #check if delete_until_special_chars has already been run
    with open(file_with_backups, "r") as checking:
        first_line = checking.readline().strip()
        if first_line.startswith("*"):
            delete_until_special_chars(file_with_backups)
        else:
            pass
        
            
    print("runnong error contin")
    sequence_letters = read_error_save_txt_file(file_path)
    sequence_numbers = convert_str_to_int_seq(sequence_letters, inverse_mapping)
    """with open(file_with_backups, "r") as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("!!" + str(progress_report_freq)):
                continue"""
    print(sequence_numbers)
    




file_with_backups = "/Users/dduren3/Desktop/nucleic_acid_replication/output/2024-02-09 22:35/Error_Contingency_run_2024-02-09 22-35.txt"

error_continuation(400, 500, 250, file_with_backups)





#ERROR CONTINUITY STILL UNDER DEV 2/16/24


#NOVEL 2 17 24

def generate_line_plot_of_longest_strand(list_w_lengths, output_folder_name, formatted_time, n_iterations):
    x_values = range(1,n_iterations +1)


    plt.plot(x_values, list_w_lengths)

    plt.xlabel("Iteration number")
    plt.ylabel("Length of longest strand")
    plt.title("Length of longest strand during each iteration")




    try:
        output_folder = "./output/" + str(output_folder_name)
        output_file_path = os.path.join(output_folder, f"_longest_strand_length_graph_{formatted_time}.jpg")
        plt.savefig(output_file_path)


    except Exception as e:
        print(e)
        print("line plot of longest strand function has failed")