"""
@author: Brandon Fujii
Final assignment for Python for Genomic Data Science @ University of Washington
"""

from sys import argv, exit
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
        print self.records
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

    def get_longest_orf(self):
        return self.max_orf
    
    def get_shortest_orf(self):
        return self.min_orf

def find_orf(record, frame = 1):
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
    for i in range(0,len(start_indices)):
        for j in range(0, len(stop_indices)):
            if start_indices[i] < stop_indices[j] and start_indices[i] > last_index:
                orf_sequence = sequence[start_indices[i]:stop_indices[j]+3]
                orfs.append(OpenReadingFrame(record.record_number, record.id, orf_sequence, start_indices[i], frame))
                last_index = stop_indices[j] + 3
                break

    return orfs

def read_fasta_file(filename):
    """ Parses a .fasta file of sequence records and returns its records """
    records = []
    for i, record in enumerate(SeqIO.parse(filename, "fasta")):
        record.record_number = i
        records.append(record)
    return records

def main():
    if len(argv) < 2:
        exit('Usage: %s filename.fasta' % argv[0])

    args = argv[1:]
    records = RecordCollection(read_fasta_file(filename = args[0]))
    orfs = records.find_all_orfs(frame = 2)
    print orfs[0], orfs[-1]

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
    
if __name__ == "__main__":
    main()