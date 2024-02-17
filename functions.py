import numpy as np
np.random.seed(0)
import seqfold
import matplotlib.pyplot as plt
import os
import time


current_time_struct = time.localtime()
formatted_time = str(time.strftime("%Y-%m-%d %H:%M", current_time_struct))


def order_of_mag(a):
    """Calculate the order of magnitude of integer elements in an array.

    The order of magnitude represents the number of nucleotides in a
    nucleic acid strand.

    Parameters
    ----------
    a : np.ndarray
        Array of integer nucleic acid sequences

    Returns
    -------
    order : np.ndarray
        Array of integers representing the length of each nucleic acid strand
    """

    order = np.vectorize(len)(np.vectorize(str)(a))

    return order

def num_nucleo(nucleotide_list):
    """Calculate the total number of nucleotides in an array.

    Can be used to check that we are not accidentally deleting any
    nucleotide when we break bonds.

    Parameters
    ----------
    nucleotide_list : np.ndarray
        Array containing all nucleic acid integer sequences

    Returns
    -------
    summ : int
        Total number of nucleotides in all strands
    """

    summ = 0
    for el in nucleotide_list:
        summ += len(str(el))
    return summ

def break_short(short, cleav_prop):
    """Randomly break bonds in short nucleic acid sequences.

    Short strands have fewer nucleotides than a given threshold.

    Parameters
    ----------
    short : np.ndarray
        Array containing all short nucleic acid integer sequences 
    cleave_prop : float
        The probability of breaking each bond

    Returns
    -------
    new_short : np.ndarray
        Array containing all broken and intact nucleic acid strands
    """

    # Calculate the number of nucleotides in each strand
    order = order_of_mag(short).astype(object)
    # Calculate the total number of bonds 
    num_bonds = np.sum(order - 1)
    # Generate an array of random numbers between 0 and 1 (length: #bonds)
    # and check whether each element is lower than a given threshold.
    # If it is lower, then we will break the bond
    cleave = np.random.random(num_bonds) < cleav_prop
    new_short = []
    i = 0 
    # Loop over all strands in the array
    for ns, seq in enumerate(short):
        part = seq 
        n_bond = 1 
        # Loop over all bonds in the strand
        for no in range(1, order[ns]):
            i += 1

            # If bond is not broken, but this is last sub-string of the
            # strand, then save
            if (not cleave[i-1]) and (no == order[ns]-1):
                new_short.append(part)
                continue

            # If bond is not broken, continue
            if not cleave[i-1]:
                n_bond += 1
                continue

            # Break bond. We represent breaking bonds by integer division
            # and modulus operation. This will separate an integer number
            # to two parts.
            part0 = part%10**n_bond
            part = part//10**n_bond
            n_bond = 1 

            # Append one part of the strand, and keep the other part
            # whose remaining bonds may be broken
            new_short.append(part0)

            # Check if the remaining part has no remaining bond, then append
            if no == order[ns]-1:
                new_short.append(part)
    
    return new_short


def structured_regions(structs):
    """Predict structured regions in a nucleic acid strand using seqfold 

    Identifies stacked base pairs, hairpins, and bulges in the sequence.

    Parameters
    ----------
    structs : list of seqfold.fold.Struct objects

    Returns
    -------
    struct_bonds : np.ndarray
        indices of bonds in structured regions
    """

    struct_bonds = []

    # Loop over all structured regions in a given strand
    for struct in structs:
        i, j = struct.ij[0][0], struct.ij[0][1] #initial and final nucleo in struct
        # If we have stacked base pairs, include indices of the bonds
        # between the stacked nucleobases
        if "STACK" in struct.desc:
            length = len(struct.desc.split(":")[1].split("/")[0])
            for bond in range(i, i+length-1):
                struct_bonds.append(bond)
            for bond in range(j-length+1, j):
                struct_bonds.append(bond) 

        # If we have hairpins or bulges, include indices of all the
        # bonds between the beginning and end of the structured region
        elif ("HAIRPIN" in struct.desc) or ("BULGE" in struct.desc):
            for bond in range(i, j):
                struct_bonds.append(bond)

            
    struct_bonds = np.array(sorted(struct_bonds), dtype=np.int64)

    return struct_bonds


def break_long(long, cleav_prop, cleav_prop_struct, mapping):
    """Randomly break bonds in long nucleic acid sequences.

    Long strands have more nucleotides than a given threshold.
    For long strands, we first identify structured regions in the strand.
    If a bond is in a structured region, we break it with probability
    `cleave_prop_struct`. If the bond is in unstructured region,
    we break it with probability `cleav_prop`, which is equal to
    the probability of breaking short strands.

    Parameters
    ----------
    long : np.ndarray
        Array containing all long nucleic acid integer sequences 
    cleave_prop : float
        The probability of breaking each bond in the unstructured region
    cleave_prop_struct : float
        The probability of breaking each bond in the structured region
    mappint : dict
        Map integer numbers (1, 2, 3, 4) with nucleobase names (A, G, C, U)

    Returns
    -------
    new_long : np.ndarray
        Array containing all broken and intact nucleic acid strands
    """

    if long.size == 0:
        return long
    
    new_longs = []

    # Convert integer sequences to base names
    long_s = convert_int_to_str_seq(long, mapping)

    # Loop over all long strands
    for seq, int_seq in zip(long_s, long):
        # Identify structured regions
        structs = seqfold.fold(seq)
        struct_bonds = structured_regions(structs)

   # Calculate the length  of the strand and the number of bonds
        order = len(seq)
        num_bonds = order - 1

        # Determine broken bonds for structured and unstructured regions
        cleave = np.random.random(num_bonds) < cleav_prop
        cleave[struct_bonds] = np.random.random(struct_bonds.size) < cleav_prop_struct

        new_long = []

        i = 0
        part = int_seq
        n_bond = 1

        # Loop over all bonds in the strand
        for no in range(1, order):
            i += 1

            # If bond is not broken, but this is last sub-string of the
            # strand, then save
            if (not cleave[i-1]) and (no == order-1):
                new_long.append(part)
                continue

            # If bond is not broken, continue
            if not cleave[i-1]:
                n_bond += 1
                continue

            # Break bond. We represent breaking bonds by integer division
            # and modulus operation. This will separate an integer number
            # to two parts.
            part0 = part%10**n_bond
            part = part//10**n_bond
            n_bond = 1

            # Append one part of the strand, and keep the other part
            # whose remaining bonds may be broken
            new_long.append(part0)

            # Check if the remaining part has no remaining bond, then append
            if (no == order-1):
                new_long.append(part)

        new_longs.extend(new_long)

    new_longs = np.array(new_longs, dtype=object)
                
    return new_longs

def convert_int_to_str_seq(seq_array, mapping):
    """Convert integer nucleic acid seqences to strings with standard names.

    Parameters
    ----------
    seq_array : np.ndarray
        Array of integer nucleic acid strands
    mapping : dict
        Map integer numbers (1, 2, 3, 4) with nucleobase names (A, G, C, U)
    """
    seq_array_s = seq_array.astype(str)

    # Replace each number with the corresponding base name
    for k, v in mapping.items():
        seq_array_s = np.char.replace(seq_array_s, k, v)

    return seq_array_s


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


"""
def generate_sequence_length_histogram(file_path, folder_name, iter_num =0):
    string_iter_num = "_" + str(iter_num)
    opened = open(file_path, "r")
    sequences = opened.read().splitlines()


    sequence_lengths = []
    for seq in sequences:
        if seq.startswith("*"):
            continue
        else:
            sequence_lengths.append(len(seq))
            

    plt.hist(sequence_lengths, bins=100)
    plt.yscale("log")
    plt.xlabel('Sequence Length (nt)')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Histogram' + str(formatted_time))

    #if file path of output is ever changed from ./output/ , below will need to be changed 
    output_folder = "./output/" + str(folder_name)
    output_file_path = os.path.join(output_folder, f"sequence_length_histogram_{formatted_time}_{string_iter_num}.jpg")
    plt.savefig(output_file_path)


    opened.close()

"""
def generate_sequence_length_histogram(file_path, folder_name, progress_report_num):

    opened = open(file_path, "r") # needs to be error contingency file
    sequences = opened.read().splitlines()
        
    sequence_lengths = []

    for i, seq in enumerate(sequences):
        if seq.startswith("!!" + str(progress_report_num)):
            break
    for seq in sequences[i+1:]:
        sequence_lengths.append(len(seq))
            
    plt.clf() # clears any data that may still be cached
    plt.hist(sequence_lengths, bins=100)
    plt.yscale("log")
    plt.xlabel('Sequence Length (nt)')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Histogram' + str(formatted_time) + "___run_" + str(progress_report_num))

    try:
        output_folder = "./output/" + str(folder_name)
        output_file_path = os.path.join(output_folder, f"{progress_report_num}_sequence_length_histogram_{formatted_time}.jpg")
        plt.savefig(output_file_path)

    except Exception as e:
        return e

    opened.close()


    
def create_blank_text_file(base_directory, folder_name, param_for_file_name): 
    #param_file name describes what the txt file being made is used for --
    #ex. error_contingency means this txt file produced is backups of each run


    #Creates a new text file with a unique name based on time function was run
    #Saves it to the location on the harddrive specified by base_directory
    
    try:
        if not os.path.exists("./output/"):
            os.mkdir("./output/")

        output_folder = os.path.join(base_directory, folder_name)
        os.mkdir(output_folder)
        file_name = f"{param_for_file_name}_run_{folder_name}.txt" #note param being used in title of file
        file_path = os.path.join(output_folder, file_name)
        with open(file_path, "w"):
            pass
        
        return file_path, output_folder
  
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return False



"""NOVEL 2 9 24"""

     
def error_safeguard_system(file_being_written_to, it_num, list_of_sequences):
    try:
        file_being_written_to.write("!!" + str(it_num) + "\n")
        for seq in list_of_sequences:
            file_being_written_to.write(seq + "\n")


    except Exception as E:
        print(E)



def generate_final_seq_file(directory, final_file_name):
    try:
        os.makedirs(directory, exist_ok=True)
        file_path = os.path.join(directory, final_file_name)
        
        with open(file_path, "w"):
            pass  
        
        print("Final text file created.")
        return file_path  
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


#NOVEL 2/17/24
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