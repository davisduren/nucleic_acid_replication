#Feed in structures overhaul


###COPY over functions.py

import numpy as np
np.random.seed(0)
import seqfold
import matplotlib.pyplot as plt
import os
import time
import sys


#last update - 3/19/24
#starting %nt involved in structured regions function

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






def structured_regions(structs, dG_critical = -6.5): 
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

    #test_list_dG_structs = []

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
  
            for bond in range(i,j+1):
                struct_bonds.append(bond)

            
        elif ("BULGE" in struct.desc):
            for bond in range(i, j):
                struct_bonds.append(bond)

            
    struct_bonds = np.array(sorted(struct_bonds), dtype=np.int64)



    #print(test_list_dG_structs)
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
        struct_bonds  = structured_regions(structs)
        #NOVEL 3/3/24- tetraloops_seen_this_iteration is just there so struct_bonds isnt a tuple and 
        #so tetraloops_seen_this_iteration doesnt break anything
        #also, this since this function is the only time structured_regions() is called
        #it must be the one to return the tetraloops_seen_this_iteration to main ipynb
        #REVISION 3/19/24 - the above, we have now decided, is not useful and 
        #this tetraloop function has since been inactivated 

   # Calculate the length  of the strand and the number of bonds
        order = len(seq)
        num_bonds = order - 1

        # Determine broken bonds for structured and unstructured regions
        cleave = np.random.random(num_bonds) < cleav_prop
        cleave[struct_bonds] = np.random.random(struct_bonds.size) < cleav_prop_struct

        new_long = []

        i = 0
        part = int_seq
        print(part)
        print("line 256")
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

    
def create_blank_text_file(base_directory, folder_name, param_for_file_name = ""): 
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
        
        #print("Final text file created.")
        return file_path  
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


    
#3/19/24 - the below has generally been decided to be less good at showing steady state
#than 
def generate_line_plot_of_longest_strand(list_w_lengths, output_folder_name, formatted_time, n_iterations):
    x_values = range(1,n_iterations +1)


    plt.clf()

    plt.yscale('linear')
    plt.ylim(0,100) # or however long we want to show
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



def percent_over_certain_length(list_of_lengths, len1, len2, len3, len4):
    """this function allows study of the effect of ratio between probability of structured and 
    unstructured regions breaking. This function takes in a list_of_lengths, which is a list
    created during every iteration that contains the lengths of each separate strand during that
    iteration. It returns the percentage of strands that are over the length dictated by
    len1, len2, etc. These percentages are then passed onto percentage_plot for plotting.
    """

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
    
    # Calculate the percentages and round to three decimal places
    percent_len_4 = round((num_len_4 / number_of_unique_strands) * 100, 3)
    percent_len_3 = round((num_len_3 / number_of_unique_strands) * 100, 3)
    percent_len_2 = round((num_len_2 / number_of_unique_strands) * 100, 3)
    percent_len_1 = round((num_len_1 / number_of_unique_strands) * 100, 3)

    #returns shortest --> longest
    return percent_len_1, percent_len_2, percent_len_3, percent_len_4
    





def percentage_plot(list_of_percentages, output_folder_name, formatted_time, len1, len2, len3, len4):
    plt.clf()
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

        #remember, list_of_percentages[0][0] correspinds to percentages relating to len1,
        #list_of_percentages[0][1] corresponds to percentages relating to len2, etc

        y_values = [sublist[i] for sublist in list_of_percentages]  # Extract the items at index i from each sublist
        plt.plot(x, y_values, label=f"Sublist Item {i+1}")

    plt.ylim(None, 100)
    plt.xlabel("Number of iterations")
    plt.ylabel("Percentage of strands over given length (%)")
    plt.title("Percentage of strands over a given length over time" + str(formatted_time))
    plt.legend(['>=' + str(len1), '>=' + str(len2), '>=' + str(len3), '>=' + str(len4)])  # Set labels for the legend

    try:
        output_folder = "./output/" + str(output_folder_name)
        output_file_path = os.path.join(output_folder, f"Percentage_strands_over_lengthX_{formatted_time}.jpg")
        plt.savefig(output_file_path, dpi=200)


    except Exception as e:
        print(e)
        print("Saving plot of percentage plot function has failed")



#BELOW NOVEL 3/19/24

def percentage_nt_involved_in_structure(list_of_sequences_being_checked, mapping):
    #list_of_sequences_being_checked is the list of the x longest seqs
    #need to find the indices of structured regions within each seq, 
    #save those numbers, and compute the percentage of each seq thag those indices constitute

    sequences = convert_int_to_str_seq(list_of_sequences_being_checked, mapping)
    #print(sequences)

    list_of_percentage_of_nt_in_structure = []
    for seq in sequences[:1]:
        try:
            structs = seqfold.fold(seq)
            struct_bonds  = structured_regions(structs)
            indices = []
            #theres a problem here - when we have structured regions that overlap,
            #like if nt at index 6 is involved in a hairpin, but also a stack, it will 
            #get counted twice -- i.e it wil go ...6 6 ...
            #need a way to just count it once
            #this is what below does

            for index in struct_bonds:
                if index not in indices:
                    indices.append(index)

            #indices contains all the number indices of the nt that are involved in structure (without repears)
            
            percentage_of_nt_in_structure = round((len(indices) / len(seq)) * 100, 3)
            list_of_percentage_of_nt_in_structure.append(percentage_of_nt_in_structure)       
    
        except Exception as e:
            print(e)
    #print(list_of_percentage_of_nt_in_structure)
    return list_of_percentage_of_nt_in_structure

    
def percentage_nt_plot(list_of_percentage_nt_involved_in_structure,output_folder_name, formatted_time):
    plt.clf()
    x_values = range(len(list_of_percentage_nt_involved_in_structure)) #should be the same as n_iterations
    y_values = list_of_percentage_nt_involved_in_structure


    plt.plot(x_values, y_values)
    plt.ylim(0,100)
    plt.xlabel("Number of iterations")
    plt.ylabel("Percentage of sequence involved w/ structure")
    plt.title("Percentage of X longest sequences that are involved in structure")

    try:
        output_folder = "./output/" + str(output_folder_name)
        output_file_path = os.path.join(output_folder, f"Nt_percentage_{formatted_time}.jpg")
        plt.savefig(output_file_path, dpi=200)


    except Exception as e:
        print(e)
        print("Saving plot of nt_percentage plot function has failed")



### END COPY OVER FUNCTIONS
        
############################################################################################

#BEGIN REPLICATION.IPYNB
        

import numpy as np
import time
import seqfold
seed_value = 0
np.random.seed(seed_value)

#LAST UPDATE: 6/24/24

#below are the variables we manipulate

init_nuc_num = 10000    # number of each base (10k A, 10k U, etc)
cleav_prop = 0.20    # chance of unstruc regions breaking during given run                  
cleav_prop_struct = 0.02    # chance of struc region breaking during given run
length_threshold = 10   # we wont check if something less than this long has struc
n_iterations = 50    # how many runs until completion
progress_report_freq  = 10   # how often code gives us a progress report / saves data to Error_Contingency file in case of error
dG_critical_value = -7.5

len1 = 4
len2 = 6    # these are the lengths we are using in the percent_over_certain_length functions
len3 = 8    #  see that function for more info
len4 = 10

list_of_percentage_nt_involved_in_structure = [] #used via percentage_nt_involved_in_structure function



#below are string versions of the variables used - this is written to our output file to help us track what settings produced the sequences
init_nuc_num_str = "*init_nuc_num = " + str(init_nuc_num) + "\n"
cleav_prop_str = "*cleav_prop = " + str(cleav_prop) + "\n"
cleav_prop_struct_str = "*cleav_prop_struct = " + str(cleav_prop_struct) + "\n"
length_threshold_str = "*length_threshold = " + str(length_threshold) + "\n"
n_iterations_str = "*n_iterations = " + str(n_iterations) + "\n"
progress_report_freq_str = "*progress_report_freq = " + str(progress_report_freq) + "\n"
seed_value_str = "*numpy random seed value = " + str(seed_value) + "\n"
settings_used = init_nuc_num_str + cleav_prop_str + cleav_prop_struct_str +  \
    length_threshold_str +  n_iterations_str + progress_report_freq_str + seed_value_str +  "\n"


#the below are used to create the file names
current_time_struct = time.localtime()
formatted_time = str(time.strftime("%Y-%m-%d %H:%M", current_time_struct))
base_directory = "./output/"
folder_name = formatted_time
blank_txt_file_paramater = "Error_Contingency" # see create_blank_text_file for description


file_path, output_directory = create_blank_text_file(base_directory, folder_name, blank_txt_file_paramater)
print(file_path)
print(output_directory)
file_to_write_to = open(file_path, "a")
file_to_write_to.write(settings_used)
final_file_path = None

list_of_percent_lens = []


###NOVEL 6/20/24
throwing_long_ones_out = True # this means we are throwing away all sequences that are longer than a set length to study the emergence of novel sequences
length_to_throw_out = 20

if throwing_long_ones_out == True:
    seqs_over_certain_length_filename = f"Sequences_over_{length_to_throw_out}_in_length.txt"
    seqs_over_certain_length_filepath = generate_final_seq_file(output_directory, seqs_over_certain_length_filename)
    with open(seqs_over_certain_length_filepath, "a") as seqs_over_certain_length_write:
        try:
            seqs_over_certain_length_write.write(settings_used)
            seqs_over_certain_length_write.flush()  # Ensure the data is flushed to the file
        except Exception as e:
            print(f"An error occurred while writing to the file: {e}")
###NOVEL 6/20/24


###NOVEL 6/27/24 feed in file
init_seq_file = "TEST_HAIRPIN_FILE.txt"
file_with_starting_seqs = open(init_seq_file, "r")
list_of_starting_seqs = []
for strand in file_with_starting_seqs:
    list_of_starting_seqs.append(strand)

nucleotide_list = [s.strip() for s in list_of_starting_seqs]

nucleotide_list = np.array(nucleotide_list)


mapping = {"1": "A", "2": "G", "3": "C", "4": "U"}


#Need to convert nucleotide_list back into list of ints for below operations

reverse_mapping = {v: k for k, v in mapping.items()}

def convert_strings_to_integers(string_nt_list):
    converted_list = []

    # Process each string in the input list
    for string in string_nt_list:
        # Convert each character in the string using the reverse mapping
        integer_string = ''.join(reverse_mapping[char] for char in string)
        # Convert the concatenated string to an integer
        integer_value = int(integer_string)
        # Add the integer value to the final list
        converted_list.append(integer_value)
    
    return converted_list


int_nt_list = convert_strings_to_integers(nucleotide_list)
for item in int_nt_list:
    print(item)                                         
### BELOW BEGINS NUCLEIC COMPUTATIONS



longest_strand_length_list = []
list_of_seqs_over_length_threshold = []

for it in range(1, n_iterations + 1):
    print(str(it) + " -- current it num")
    
    strand_length_tracking_list = []

    # Phase 1: Pair nucleotide strands
    # Shuffle list of nucleotide strands
    np.random.shuffle(nucleotide_list)

    # Take the first half of the strands and calculate their lengths
    size = nucleotide_list.size
    # Use dtype=object to keep arbitrary integer prevision
    order = order_of_mag(nucleotide_list[:size//2]).astype(object)


    # Phase 2: Determine folded structures and randomly break bonds in long strands
    order = order_of_mag(nucleotide_list).astype(object)
    mono = nucleotide_list[order == 1]
    short = nucleotide_list[np.logical_and(order > 1, order < length_threshold)]
    long = nucleotide_list[order >= length_threshold]
    
    

    long = break_long(long, cleav_prop, cleav_prop_struct, mapping)

        # Phase 3: Randomly break bonds in short strands
    short = break_short(short, cleav_prop)

    nucleotide_list = np.concatenate((mono, short, long))
    




    strands_in_decreasing_order = np.sort(nucleotide_list)[::-1]
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

            


        
    #this stuff below makes a list of the percentages over certain lengths 
    #it appends the list of numbers (which are percentages) to another list, list_of_percent_lens
    percent_len1, percent_len2, percent_len3, percent_len4 = percent_over_certain_length(strand_length_tracking_list, len1,len2, len3, len4)
    sublist_percent_lens_during_each_it = [percent_len1, percent_len2, percent_len3, percent_len4]
    list_of_percent_lens.append(sublist_percent_lens_during_each_it)




    longest_strand = strands_in_decreasing_order[0]
    longest_strand_str = str(longest_strand)
    longest_strand_length_list.append(len(longest_strand_str))


    #below is 3/19/24
    summary_for_percent_nt = np.sort(nucleotide_list)[::-1][:10] #gives the 10 longest sequences
    
    list_of_percentage_nt_involved_in_structure.append(percentage_nt_involved_in_structure(summary_for_percent_nt, mapping))


    if it%n_iterations == 0: 
        #means code is done, generate final txt file with end results
        final_file_name = f"FINAL_run_{folder_name}.txt"
        final_file_path = generate_final_seq_file(output_directory, final_file_name)
        final_file_to_write_to = open(final_file_path, "a")
        final_file_to_write_to.write(settings_used)
        final_file_to_write_to.write("!!" + str(n_iterations) + "\n")

        summary = np.sort(nucleotide_list)[::-1] #to get all lengths, remove [::-1][-:10]
        summary_s = convert_int_to_str_seq(summary, mapping)


        for s in summary_s:
            final_file_to_write_to.write(s + "\n")


         #throwing_out_long_seqs:
            
        if throwing_long_ones_out == True:
            list_of_seqs_over_length_threshold = np.array(list_of_seqs_over_length_threshold)
            str_list_of_long_seqs = convert_int_to_str_seq(list_of_seqs_over_length_threshold, mapping)
            with open(seqs_over_certain_length_filepath, "a") as seqs_over_certain_length_write:
                for s in str_list_of_long_seqs:
                    print(s)
                    seqs_over_certain_length_write.write(s + "\n")
        


    elif it%progress_report_freq == 0:
        print("\nRESULTS\n")

        summary = np.sort(nucleotide_list)[::-1] # [::-1] = longest --> shortest
        print(summary)
        summary_s = convert_int_to_str_seq(summary, mapping)

        list_of_seqs = []
        for s in summary_s:
            list_of_seqs.append(s)
        error_safeguard_system(file_to_write_to, it, list_of_seqs)

        generate_sequence_length_histogram(file_path, folder_name, it)





generate_sequence_length_histogram(final_file_path, folder_name, n_iterations)
#generate_line_plot_of_longest_strand(longest_strand_length_list, folder_name, formatted_time, n_iterations)
#print(list_of_percent_lens)
percentage_plot(list_of_percent_lens, folder_name, formatted_time, len1, len2, len3, len4)
#percentage_nt_plot(list_of_percentage_nt_involved_in_structure, folder_name, formatted_time)

