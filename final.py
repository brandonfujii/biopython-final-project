"""
@author: Brandon Fujii
Final assignment for Python for Genomic Data Science @ University of Washington
"""

from sys import argv, exit
from Bio import SeqIO

class OpenReadingFrame(object):
    def __init__(self, record_number, record_id, sequence, start, frame = 1):
        self.record_number = record_number
        self.record_id = record_id
        self.sequence = sequence
        self.start = start
        self.frame = frame

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return "<OpenReadingFrame record_id = {}, record_number = {}, sequence_length = {}, start_position = {}, frame = {}>".format(self.record_id, self.record_number, self.__len__(), self.start, self.frame)

class RecordCollection(object):
    def __init__(self, records):
        self.records = dict([(record.id, record) for record in records])
        self.max_sequences = self.get_longest_sequence_records()
        self.min_sequences = self.get_shortest_sequence_records()
        self.orfs = []
    
    def __len__(self):
        return len(self.records)

    def get_record_by_id(self, id):
        try:
            record = self.records[id]
            return record
        except KeyError as e:
            return None

    def get_longest_sequence_records(self):
        max_sequence_length_records = []

        for id, record in self.records.items():
            sequence_length = len(record.seq)
            if not max_sequence_length_records or len(max_sequence_length_records[0].seq) == sequence_length:
                max_sequence_length_records.append(record)
            elif sequence_length > len(max_sequence_length_records[0].seq):
                max_sequence_length_records = [record]

        self.max_sequences = max_sequence_length_records
        return max_sequence_length_records 

    def get_shortest_sequence_records(self):
        min_sequence_length_records = []

        for id, record in self.records.items():
            sequence_length = len(record.seq)
            if not min_sequence_length_records or len(min_sequence_length_records[0].seq) == sequence_length:
                min_sequence_length_records.append(record)
            elif sequence_length < len(min_sequence_length_records[0].seq):
                min_sequence_length_records = [record]

        self.min_sequences = min_sequence_length_records
        return min_sequence_length_records 

    def find_all_orfs(self, frame = 1):
        orfs = []
        for id, record in self.records.items():
            if record.seq:
                record_orfs = find_orf(record, frame)
                orfs.extend(record_orfs)

        self.orfs = orfs
        return sorted(orfs, key = lambda x: len(x), reverse=True)
    
    def find_all_repeats(self, repeat_length):
        repeats = {}
        for id, record in self.records.items():
            if record.seq:
                repeats = find_repeats(record.seq, repeat_length, repeats)
        
        return repeats

def find_orf(record, frame = 1):
    """
    Finds a list of open reading frames, which is a sub-sequence that begins with a start codon ("ATG"), and ends with\
    an end codon ("TAA", "TGA", "TAG"), based on a forward reading frame
    """
    start_indices = []
    stop_indices = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TGA", "TAG"]
    sequence = record.seq

    for i in range(frame - 1, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon == start_codon:
            start_indices.append(i)

        if codon in stop_codons:
            stop_indices.append(i)

    orfs = []
    last_index = 0
    for start_index in start_indices:
        for stop_index in stop_indices:
            if start_index < stop_index and start_index >= last_index:
                orf_sequence = sequence[start_index:stop_index + 3]
                orfs.append(OpenReadingFrame(record.record_number, record.id, orf_sequence, start_index, frame))
                last_index = stop_index + 3
                break

    return orfs

def read_fasta_file(filename):
    """ Parses a .fasta file of sequence records and returns its records """
    records = []
    for i, record in enumerate(SeqIO.parse(filename, "fasta")):
        record.record_number = i
        records.append(record)
    return records

def find_repeats(sequence, repeat_length, repeats_dict = {}):
    """
    Returns a dictionary of repeats (repeat_sequence -> num_occurrences) of a specified length
    in a given sequence
    """
    for i in range(len(sequence)):   
        read = str(sequence[i:i + repeat_length])
        if read and len(read) == repeat_length:
            if read not in repeats_dict:
                repeats_dict[read] = 1
            else:
                repeats_dict[read] += 1
    return repeats_dict

def find_max_in_dict(dictionary):
    return max(zip(dictionary.values(), dictionary.keys())) if len(dictionary) else None

def main():
    if len(argv) < 2:
        exit('Usage: %s filename.fasta' % argv[0])

    args = argv[1:]
    records = RecordCollection(read_fasta_file(filename = args[0]))

    print "============ Question 1 ============"
    """
    Q1: How many records are in the multi-FASTA file?
    """
    print len(records)

    print "============ Question 2 ============"
    """
    Q2: What is the length of the longest sequence in the file?
    """
    longest_sequences = records.get_longest_sequence_records()
    print len(longest_sequences[0]) if longest_sequences else 0

    print "============ Question 3 ============"
    """
    Q3: What is the length of the shortest sequence in the file?
    """
    shortest_sequences = records.get_shortest_sequence_records()
    print len(shortest_sequences[0]) if shortest_sequences else 0

    print "============ Question 4 ============"
    """
    Q4: What is the length of the longest ORF appearing in reading frame 2 of any of the sequences?
    """
    q4_orfs = records.find_all_orfs(frame = 2)
    print len(q4_orfs[0])

    print "============ Question 5 ============"
    """
    Q5: What is the starting position of the longest ORF in reading frame 3 in any of the sequences?
    The position should indicate the character number where the ORF begins.
    """
    q5_orfs = records.find_all_orfs(frame = 3)
    print q5_orfs[0].start if q5_orfs else None

    print "============ Question 6 ============"
    """
    Q6: What is the length of the longest ORF appearing in any sequence and in any forward reading frame?
    """
    max_orf_length = 0
    for frame in range(1, 4):
        q6_orfs = records.find_all_orfs(frame = frame)
        if len(q6_orfs[0]) > max_orf_length:
            max_orf_length = len(q6_orfs[0])

    print max_orf_length

    print "============ Question 7 ============"
    """
    Q7: What is the length of the longest forward ORF that appears in the sequence with the identifier gi|142022655|gb|EQ086233.1|16?
    """
    max_orf_length = 0
    for frame in range(1, 4):
        q7_orfs = find_orf(records.get_record_by_id("gi|142022655|gb|EQ086233.1|16"), frame = frame)
        max_length = max(len(orf) for orf in q7_orfs)

        if max_length > max_orf_length:
            max_orf_length = max_length

    print max_orf_length

    print "============ Question 8 ============"
    """
    Q8: Find the most frequently occurring repeat of length 6 in all sequences. How many times does it occur in all?
    """
    repeats = records.find_all_repeats(repeat_length = 6)
    num_occurrences, most_occurring_repeat = find_max_in_dict(repeats)
    print most_occurring_repeat, "occurs", num_occurrences, "times"

    print "============ Question 9 ============"
    """
    Q9: Find all repeats of length 12 in the input file. Let's use Max to specify the number of copies
    of the most frequent repeat of length 12. How many different 12-base sequences occur Max times?
    """
    repeats = records.find_all_repeats(repeat_length = 12)
    num_occurrences, most_occurring_repeat = find_max_in_dict(repeats)
    max_repeats = []
    for key, value in repeats.items():
        if value == num_occurrences:
            max_repeats.append(key)
    print len(max_repeats)

    print "============ Question 10 ============"
    """
    Q10: Which one of the following repeats of length 7 has a maximum number of occurrences?
    """
    repeats = records.find_all_repeats(repeat_length = 7)
    num_occurrences, most_occurring_repeat = find_max_in_dict(repeats)
    print most_occurring_repeat

if __name__ == "__main__":
    main()