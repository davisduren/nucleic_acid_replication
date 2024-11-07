import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seqfold
import os
from functions import *
import time
import sys
from seqfold import dg, dg_cache, fold, Cache, Struct, dot_bracket
from collections import defaultdict







"""def process_cleaved_sequences_with_labels(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    sequence_fragments = {}
    current_label = None
    original_sequences = []

    for line in lines:
        line = line.strip()
        if line.startswith("!!"):  # New sequence section
            current_label = int(line[2:])  # Extract sequence number
            sequence_fragments[current_label] = []  # Initialize the list for this sequence
        elif line and not line.startswith("*"):  # Process valid sequence lines
            sequence_fragments[current_label].append(line)

    return sequence_fragments


def count_cleavage_frequencies(sequence_fragments, original_sequences):
    cleavage_frequencies = [0] * (len(original_sequences[0]) - 1)  # Create a list for positions
    
    # Track cleaved positions
    for seq_label, fragments in sequence_fragments.items():
        if seq_label <= len(original_sequences):  # Ensure seq_label is within bounds
            original_seq = original_sequences[seq_label - 1]  # Get the corresponding original sequence
            for fragment in fragments:
                # Compare the original sequence with the fragment
                for index in range(len(original_seq) - 1):  # Check each bond position
                    if index < len(fragment) and original_seq[index] != fragment[index]:  # If there's a change, cleavage occurred
                        cleavage_frequencies[index] += 1

    return cleavage_frequencies


# Example usage
original_sequences = [
    "AAAAAAGGGGGGUUUUCCCCCC",  # Original Sequence 1
    "AAAAAAGGGGGGUUUUCCCCCC",  # Original Sequence 2
    "AAAAAAGGGGGGUUUUCCCCCC",  # Original Sequence 3
]

# Path to your file containing the cleaved sequences
file_path = "/Users/dduren3/Desktop/nucleic_acid_replication/output/2024-09-30 17:21/Error_Contingency_run_2024-09-30 17:21.txt"  # Update with your file path

# Process sequences
sequence_fragments = process_cleaved_sequences_with_labels(file_path)
print("Processed Sequence Fragments:")
for label, fragments in sequence_fragments.items():
    print(f"Label {label}: {fragments}")

# Count cleavage frequencies
cleavage_frequencies = count_cleavage_frequencies(sequence_fragments, original_sequences)
print("Cleavage Frequencies:")
for position, count in enumerate(cleavage_frequencies):
    print(f"Position {position}: Cleaved {count} times.")
"""
#10/16/24 update




def calculate_cleavage_probabilities(file_path, original_sequence, num_of_orig_seq):
    # Read the output file containing sequences after the operation
    with open(file_path, 'r') as f:
        processed_sequences = f.readlines()

    processed_sequences = [
        seq.strip() for seq in processed_sequences 
        if seq.strip() and not seq.startswith('*') and not seq.startswith('!!')
    ]

    total_sequences = num_of_orig_seq 

    seq_length = len(original_sequence)

    # Initialize a list to track the number of times each bond was cleaved
    cleavage_counts = [0] * (seq_length - 1)  





    total_bonds_counted = [0] * (seq_length - 1)

    # Process each sequence and analyze the cleavage points
    i = 0
    while i < len(processed_sequences):
        fragments = []  # Reset the fragments for each new sequence group
        if processed_sequences[i] == original_sequence:
            # If the sequence is intact, count all bonds as unbroken
            print(f"Found intact sequence: {processed_sequences[i]}")

            for bond_index in range(seq_length - 1):
                total_bonds_counted[bond_index] += 1
            i += 1
            continue

        # Collect fragments until a new sequence group starts
        while i < len(processed_sequences) and processed_sequences[i] != original_sequence:
            fragments.append(processed_sequences[i])
            i += 1

        # If we have multiple fragments, analyze the cleavage point
        if fragments:
            print(f"Fragments found: {fragments}")
            for frag in fragments:
                frag = frag.strip()
                frag_length = len(frag)

                # Find the matching index in the original sequence starting from start_idx
                for start in range(len(original_sequence) - len(frag) + 1):
                    if original_sequence[start:start + len(frag)] == frag:
                        match_index = start
                        break

                if match_index != -1:
                    # Correct bond index where the fragment ends
                    if match_index + frag_length < seq_length:
                        bond_index = match_index + frag_length - 1
                        cleavage_counts[bond_index] += 1

                    # Increment total bonds counted for each bond in this fragment
                    for bond_index in range(match_index, match_index + frag_length - 1):
                        total_bonds_counted[bond_index] += 1

                    # Update start_idx to be just past the current match
                else:
                    print(f"Fragment {frag} does NOT match original sequence!")



    # Calculate cleavage probabilities for each bond
    cleavage_probabilities = [
        cleavage_counts[i] / total_sequences if total_bonds_counted[i] > 0 else 0.0
        for i in range(seq_length - 1)
    ]
    
    return cleavage_probabilities





original_sequence = "AAAAAAGGGGGGUUUUCCCCCC"
file_path = "/Users/dduren3/Desktop/nucleic_acid_replication/output/5000 1 it/FINAL_run_2024-11-07 14:41.txt"
num_of_orig_seq = 5000 # this is how many of original_sequence were used in the original simulation
cleavage_probabilities = calculate_cleavage_probabilities(file_path, original_sequence, num_of_orig_seq)
print(cleavage_probabilities)




def calculate_percentage_of_whole_seqs(file_path, original_sequence, num_of_orig_seq, num_struct, num_unstruct):


    with open(file_path, 'r') as f:
        processed_sequences = f.readlines()

    star_lines = [seq.strip() for seq in processed_sequences if seq.startswith('*')]

    processed_sequences = [
        seq.strip() for seq in processed_sequences 
        if seq.strip() and not seq.startswith('*') and not seq.startswith('!!')
    ]

    cleavage_props = {}
    for line in star_lines:
        if line.startswith('*cleav_prop') or line.startswith('*cleav_prop_struct'):
            key, value = line.split('=')
            cleavage_props[key.strip()] = float(value.strip())

    cleav_prop = cleavage_props.get('*cleav_prop')
    cleav_prop_struct = cleavage_props.get('*cleav_prop_struct')

    theoretical_percentage_unbroken_seqs = round(((((1-cleav_prop)**num_unstruct) * ((1-cleav_prop_struct)**num_struct))) * 100, 2)




    num_unbroken = 0
    for seq in processed_sequences:
        if len(seq) == len(original_sequence):
            num_unbroken +=1

    exp_percentage_unbroken_seqs = round((num_unbroken / num_of_orig_seq) * 100, 2)

    return f"theoretical percentage unbroken: {theoretical_percentage_unbroken_seqs}, " \
       f"experimental percentage unbroken: {exp_percentage_unbroken_seqs}"






original_sequence = "AAAAAAGGGGGGUUUUCCCCCC"
file_path = "/Users/dduren3/Desktop/nucleic_acid_replication/output/500 1 it/FINAL_run_2024-10-27 12:46.txt"
num_of_orig_seq = 500 
num_struct = 15
num_unstruct = 6

print(calculate_percentage_of_whole_seqs(file_path, original_sequence, num_of_orig_seq, num_struct, num_unstruct))






#10/8/24 update:
#below i think is working - gives color map of all indices - it just needs to be given a list 
    #containing cleavage_probs, which is what above is trying to do

def bond_cleavage_color_map(cleavage_probs):
    colors = []
    for prob in cleavage_probs:
        # Scale cleavage probability to red-blue gradient
        red_intensity = prob / 100
        blue_intensity = 1 - red_intensity
        colors.append((red_intensity, 0, blue_intensity))  
    return colors

# Test data
sequence = "AAAAAAGGGGGGUUUUCCCCCC"  
cleavage_probs = cleavage_probabilities


colors = bond_cleavage_color_map(cleavage_probs)

fig, ax = plt.subplots(figsize=(12, 3))
# Create ticks only for the bonds (sequence length - 1)
ax.set_xticks(np.arange(len(sequence) - 1) + 0.5)
# Label each bond with the residues it connects and the cleavage percentage
ax.set_xticklabels([f'{sequence[i]} - {sequence[i+1]} ({cleavage_probs[i]}%)' for i in range(len(sequence)-1)])
ax.set_yticks([])  

# Draw colored rectangles between each pair of residues
for i in range(len(sequence) - 1):
    rect = plt.Rectangle((i, 0), 1, 1, facecolor=colors[i])
    ax.add_patch(rect)

red_patch = mpatches.Patch(color='red', label='100% Cleavage')
blue_patch = mpatches.Patch(color='blue', label='0% Cleavage')
ax.legend(handles=[red_patch, blue_patch], loc='upper right')

plt.xlim(0, len(sequence) - 1)
plt.ylim(0, 1)
plt.title('Bond Cleavage Probability Visualization with Percentages')
plt.show()





"""
print(dg('AAAAAAGGGGGGUUUUCCCCCC'))

strand = 'AAAAAAGGGGGGUUUUCCCCCC'

fold = seqfold.fold(strand)


struct_bonds = []
for item in fold:
    i, j = item.ij[0][0], item.ij[0][1]
    print(item)
    if "STACK" in item.desc:
        
        length = len(item.desc.split(":")[1].split("/")[0])
        for bond in range(i, i+length-1):
            struct_bonds.append(bond)
        for bond in range(j-length+1, j):
            print(j-length+1)
            struct_bonds.append(bond) 
        

        struct_bonds.append(i)
        struct_bonds.append(j)
    if ("HAIRPIN" in item.desc):

        for bond in range(i,j+1): #should this be j+1? or just j??
            struct_bonds.append(bond)

struct_bonds.sort()
print(struct_bonds)

"""



sys.exit()


print(len(str(124312431423213420)))

init_seq_file = "TEST_HAIRPIN_FILE copy.txt"
file_with_starting_seqs = open(init_seq_file, "r")
list_of_starting_seqs = []
for strand in file_with_starting_seqs:
    list_of_starting_seqs.append(strand)

feeding_in_knowns = True

int_nt_list = [s.strip() for s in list_of_starting_seqs]

int_nt_list = np.array(int_nt_list)

mapping = {"1": "A", "2": "G", "3": "C", "4": "U"}


#Need to convert int_nt_list back into list of ints for below operations

reverse_mapping = {v: k for k, v in mapping.items()}

def convert_strings_to_integers(string_nt_list):
    converted_list = np.array([], dtype = np.float128) 

    # Process each string in the input list
    for string in string_nt_list:
        # Convert each character in the string using the reverse mapping
        integer_string = ''.join(reverse_mapping[char] for char in string)
        integer_value = int(integer_string)
        converted_list = np.append(converted_list, integer_value)
    
    return converted_list


"""int_nt_list = np.array(convert_strings_to_integers(int_nt_list), dtype=np.float128)
convert = np.vectorize(np.format_float_positional)

print(np.format_float_positional(int_nt_list))"""

int_nt_list = np.array(convert_strings_to_integers(int_nt_list), dtype=np.float128)

formatted_list = np.vectorize(lambda x: np.format_float_positional(x, precision=0, trim='-'))(int_nt_list)

print(formatted_list)
    
    


### BELOW BEGINS NUCLEIC COMPUTATIONS










##
# #6/29/24 backup
"""    strands_in_decreasing_order = np.sort(nucleotide_list)[::-1]
    for strand in strands_in_decreasing_order:
        strand_length_tracking_list.append(len(str(strand)))
        #below is NOVEL 6/20/24
        #i already had this strands in decreasing order, may as well recycle it for the "if strand is >= x length, throw it out of the pool" test

            
        if throwing_long_ones_out == True:

            if len(str(strand)) >= length_to_throw_out: # means we want to throw that strand out
                strand_numpyarray = np.array(strand)
                strand_str = str(convert_int_to_str_seq(strand_numpyarray, mapping))
                dg_value = seqfold.dg(strand_str)
                if seqfold.dg(strand_str) <= dG_critical_value: # means it is more negative than our dG_crit
                    print(strand_str + str(dg_value))
                    indices_to_remove = np.where(nucleotide_list == strand)[0]
                    list_of_seqs_over_length_threshold.append(strand) # see "if it%n_iterations ==0:" line
                    nucleotide_list = np.delete(nucleotide_list , indices_to_remove) 
                    stringified_nt_list = map(str, nucleotide_list)
                    concatenated_digits = ''.join(stringified_nt_list)
                    total_nt = len(concatenated_digits)
                    print(total_nt) # without any deletions, should be 4* whatever init_nuc_num is. If not, means this code is working

            """


sys.exit()
###6/21/24 
"""How can we quantitate the number of base pairs of a folded 20 mer or similar length sequence?"""



def structured_regions(structs, dG_critical = -1.5): 
    #if I change dG_critical to a negative value, it deletes internal loops from structs. Why? No idea...
    ###NOVEL 6/22/24 - shifting this over to dG calculation


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
        string_struct_split = str(struct).split()
        if len(string_struct_split) >=3:
            dG_value = float(string_struct_split[2])
        else:
            raise ValueError
        
        if dG_value >= dG_critical: #means that we have an unstable structure, remove from list of potential "structs"
            continue


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
                
        elif ("HAIRPIN" in struct.desc):
            i +=1
            j +=1 ### !!! CRITICAL FIX !!! 6/22/24    i think        

            for bond in range(i,j+1):
                struct_bonds.append(bond)

            
        elif ("BULGE" in struct.desc):
            for bond in range(i, j):
                struct_bonds.append(bond)

            
    struct_bonds = np.array(sorted(struct_bonds), dtype=np.int64)


    return struct_bonds


#

rna_sequence = "AGCACAAGUCUUCGCAAUGGUUUCUCUUG" 


structs = seqfold.fold(rna_sequence)
struct_bonds = structured_regions(structs)














sys.exit()


\
"""below is modified structured_regions based off of dG for each structure
def structured_regions(structs, dG_critical = -6.5): 
    #if I change dG_critical to a negative value, it deletes internal loops from structs. Why? No idea...
    ###NOVEL 6/22/24 - shifting this over to dG calculation


    Predict structured regions in a nucleic acid strand using seqfold 

    Identifies stacked base pairs, hairpins, and bulges in the sequence.

    Parameters
    ----------
    structs : list of seqfold.fold.Struct objects

    Returns
    -------
    struct_bonds : np.ndarray
        indices of bonds in structured regions
    

    struct_bonds = []

    test_list_dG_structs = []

    # Loop over all structured regions in a given strand
    for struct in structs:
        string_struct_split = str(struct).split()
        if len(string_struct_split) >=3:
            dG_value = float(string_struct_split[2])
        else:
            raise ValueError
        
        if dG_value >= dG_critical: #means that we have an unstable structure, skip it 

            continue
        if dG_value <= dG_critical:
            print(dG_value)



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
                
        elif ("HAIRPIN" in struct.desc):
            i +=1
            j +=1 ### !!! CRITICAL FIX (i think) !!! 6/22/24           

            for bond in range(i,j+1):
                struct_bonds.append(bond)

            
        elif ("BULGE" in struct.desc):
            for bond in range(i, j):
                struct_bonds.append(bond)

            
    struct_bonds = np.array(sorted(struct_bonds), dtype=np.int64)



    print(test_list_dG_structs)
    return struct_bonds"""


##

"""Below is untouched structured_regions for contingency

def structured_regions(structs):
    #below start quotations
    Predict structured regions in a nucleic acid strand using seqfold 

    Identifies stacked base pairs, hairpins, and bulges in the sequence.

    Parameters
    ----------
    structs : list of seqfold.fold.Struct objects

    Returns
    -------
    struct_bonds : np.ndarray
        indices of bonds in structured regions
    #end quotations
    

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

"""



"""NOVEL 3/1/24
Make function that will track seen tetraloops -
attempting to modify structured_regions function to do this 
"""

#below returns just a list with all possible tetraloop combos
#we use this to check later on if our tetraloop has been found yet at any time

def generate_all_combinations():
    combinations = []
    nucleotides = ['A', 'U', 'G', 'C']
    for nucleotide1 in nucleotides:
        for nucleotide2 in nucleotides:
            for nucleotide3 in nucleotides:
                for nucleotide4 in nucleotides:
                    combination = nucleotide1 + nucleotide2 + nucleotide3 + nucleotide4
                    combinations.append(combination)
    return combinations



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

    all_possible_tetraloops = generate_all_combinations()
    tetraloops = []

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
        elif ("HAIRPIN" in struct.desc):
            if j-i ==5: # means tetraloop present
                for bond in range(i,j):
                    tetraloops.append(bond)
                    print(tetraloops)
                    print("tetealoop shit working fine")

            else:
                for bond in range(i,j):
                    struct_bonds.append(bond)



            pass
            
        elif ("BULGE" in struct.desc):
            for bond in range(i, j):
                struct_bonds.append(bond)

            
    struct_bonds = np.array(sorted(struct_bonds), dtype=np.int64)

    return struct_bonds





"""NOVEL 2/23/24
Make percentage over certain length tracking function"""
def percent_over_certain_length(list_of_lengths, len1, len2, len3, len4):
    #list of lengths = list with lengths of all strands - this list is reset during each iteration
    #file_path = file that this data with the percentages will be saved to during each iteration
    #folder name = where we will save the graph and file with data to
    #len1 = shortest length we will check (ex. 10) - these are ints
    #len2 = 2nd shortest length we are interested in (ex. 20)
    #etc


    number_of_unique_strands = 0
    
    num_len_1 = 0
    num_len_2 = 0
    num_len_3 = 0
    num_len_4 = 0


    for strand in list_of_lengths:
        number_of_unique_strands +=1
        if strand >=len4:
            num_len_4 +=1
        elif strand >=len3:
            num_len_3 +=1
        elif strand >=len2:
            num_len_2 +=1
        elif strand >=len1:
            num_len_1 +=1
    
    # percentages
    # Calculate the percentages and round to three decimal places
    percent_len_4 = round((num_len_4 / number_of_unique_strands) * 100, 3)
    percent_len_3 = round((num_len_3 / number_of_unique_strands) * 100, 3)
    percent_len_2 = round((num_len_2 / number_of_unique_strands) * 100, 3)
    percent_len_1 = round((num_len_1 / number_of_unique_strands) * 100, 3)


    return percent_len_1, percent_len_2, percent_len_3, percent_len_4
    #returns shortest --> longest

    
    
def percentage_plot(list_of_percentages, output_folder_name, formatted_time):

    #list of percentages is a list of lists -- the total number of items in the main list 
    #is the same as n_iterations completed
    #each sublist is the percentage of strands over their given length at that iteration
    #first item in each sublist is len1, 2nd is len2, etc
    #formatted time = the output directory (since the name of output directory is 
    #simply the formatted time variable from main file)

    x = range(len(list_of_percentages))  # x axis will be number of iterations completed
    y_values = [[] for _ in range(len(list_of_percentages[0]))]  # makes empty lists for each y-value
    
    # Extract y-values from the sublists - each sublist is the list of percentages from any given iteration
    for i in range(len(list_of_percentages[0])):
        y_values = [sublist[i] for sublist in list_of_percentages]  # Extract the items at index i from each sublist
        plt.plot(x, y_values, label=f"Sublist Item {i+1}")
    
    plt.xlabel("Number of iterations")
    plt.ylabel("Percentage of strands over given length")
    plt.title("Percentage of strands over a given length over time")
    plt.legend()  # Add legend to show which line corresponds to which y-value
    plt.show()

    try:
        output_folder = "./output/" + str(output_folder_name)
        output_file_path = os.path.join(output_folder, f"_percentage_strands_over_lengthX_{formatted_time}.jpg")
        plt.savefig(output_file_path)


    except Exception as e:
        print(e)
        print("Saving plot of percentage plot function has failed")

    







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
def average_sequence_length_and_longest_seq(file_path):
    opened = open(file_path, "r")
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
def generate_sequence_length_histogram(file_path):
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