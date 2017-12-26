from sys import argv, exit
from Bio import SeqIO

class Sequence(object):
    "Simplified DNA sequence formatting object" 
    def __init__(self, record, id, seq):
        self.record = record
        self.id = id
        self.data = seq
        self.len = len(seq) if seq else 0
        self.orf = None
        self.orfs = []
        self.orflen = 0
        self.start_position = None
        self.repeats = {}
        self.max_repeat = 0
    
    def __repr__(self):
        if self.orf:
            return "<Sequence id = %s, len = %d, orflen = %d, orf_start_position = %d>" % (self.id, self.len, self.orflen, self.start_position)
        else:
            return "<Sequence id = %s, len = %d>" % (self.id, self.len)
    
    def summarize(self):
        if not self.orf:
            self.find_orf()

        if not self.max_repeat:
            self.find_max_repeat()

        print """\
        =================== RECORD %d ===================
        id..................%s,
        sequence length.....%d,
        open reading frame..%s,
        ORF length..........%d,
        Repeats.............%r
        Max repeat..........(%s, %d)
        """ % \
        (self.record,
         self.id, 
         self.len,
         self.orf or "N/A",
         self.orflen,
         self.repeats,
         self.max_repeat[1],
         self.max_repeat[0])

    def find_longest_orf(self):
        if not len(self.orfs):
            return None

        maxlen = []
        for orf in self.orfs:
            if not maxlen or orf['len'] == maxlen[0]['len']:
                maxlen.append(orf)
            elif maxlen[0]['len'] < orf['len']:
                maxlen = [orf]
            else:
                pass

        return maxlen 

    def find_orf(self, frame = 1):
        "Given a sequence object and a frame, finds open reading frame that begins with a start codon\
        and ends with one of three stop codons. Returns (starting_position, ORF sequence)"
        seq = str(self.data)
        seq = "ATGAAATAG"
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]

        is_orf = False
        orf = {}
        orf_seq = []
        orfs = []

        for i in range(frame - 1, len(seq) + 1, 3):
            read = seq[i:i+3]
            print read
            if read == start_codon:
                orf["start"] = i
                is_orf = True

            if is_orf:
                orf_seq.append(read)
            
            if read in stop_codons and is_orf:
                is_orf = False
                orf["sequence"] = "".join(orf_seq)
                orf["len"] = len(orf_seq)
                orfs.append(orf)
                orf_seq = []

        self.orfs = orfs
        print self.orfs
        return orfs

    def find_max_repeat(self, repeat_length=3):
        "Finds the maximum repetition of a sequence of a given length"
        if not len(self.repeats):
            self.find_repeats(repeat_length)

        self.max_repeat = max(zip(self.repeats.values(), self.repeats.keys())) if len(self.repeats) else None
        return self.max_repeat

    def find_repeats(self, repeat_length=3):
        self.repeats = {}
        seq = str(self.data)
        for i in range(len(seq)):   
            read = seq[i:i + repeat_length]
            if read and len(read) == repeat_length:
                if read not in self.repeats:
                    self.repeats[read] = 1
                else:
                    self.repeats[read] += 1
        return self.repeats

def find_extreme_length(seq, length_list, max = True):
    id = seq.id
    length = seq.len

    if not length_list or length_list[0].len == length:
        length_list.append(seq)
    elif (max and length_list[0].len < length) or (not max and length_list[0].len > length):
        length_list = [seq]
    else:
        pass

    return length_list

def find_longest_orf(seq, orf_list):
    id = seq.id
    orf = seq.orfs

    if not orfs:
        return orf_list

    if not orf_list or orf_list[0]['len'] == orfs[0]['len']:
        orf_list.append(seq)
    elif orf_list[0]['len'] < orfs[0]['len']:
        orf_list = orfs
    else:
        pass

    return orfs_list

def main():
    if len(argv) < 2:
        exit('Usage: %s filename.fasta' % argv[0])

    args = argv[1:]
    max_length_seq = []
    min_length_seq = []
    max_length_orf = []
    min_length_orf = []
    num_records = 0
    max_repeat = None

    for i, seq_record in enumerate(SeqIO.parse(args[0], "fasta")):
        num_records += 1
        seq = Sequence(i, seq_record.id, seq_record.seq)
        seq.find_orf()
        seq.find_longest_orf
        # seq.find_max_repeat(6)
        seq.summarize()
        break

        max_length_seq = find_extreme_length(seq, max_length_seq, max = True)
        min_length_seq = find_extreme_length(seq, min_length_seq, max = False)


        if not max_repeat or seq.max_repeat[0] > max_repeat[0]:
            max_repeat = seq.max_repeat

    print "\n\n=================== SUMMARY ==================="
    print "Number of records:", num_records
    print "Longest records:", max_length_seq
    print "Shortest records:", min_length_seq
    print "Longest ORFs:", max_length_orf
    print "Shortest ORFs:", min_length_orf
    print "Max repeat:", max_repeat




if __name__ == "__main__":
    main()