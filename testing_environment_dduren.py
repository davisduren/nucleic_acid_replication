"""def show_most_developed_sequences(filename):

    opened = open(filename, "r")
    listlines = opened.readlines()
    #print(listlines)
    list_of_most_developed = []
    for seq in listlines[-10:]:
        x = seq.strip("\n")
        list_of_most_developed.append(x)
    opened.close()

    appending = open(filename, "a")
    appending.write("Most developed sequences: " + "\n" + list_of_most_developed)
    appending.close()"""



#print(show_most_developed_sequences("test_file_1.txt"))


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


print(average_sequence_length_and_longest_seq("test_file_1.txt"))

print(average_sequence_length("test_file_1.txt"))

