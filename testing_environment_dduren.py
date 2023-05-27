import matplotlib.pyplot as plt





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
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Histogram')
    plt.savefig('sequence_length_histogram.jpg')  # Save the histogram as a .jpg file
    plt.show()

    opened.close()
# Usage example
generate_sequence_length_histogram('test_file_1.txt')

