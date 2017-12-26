import unittest
from random import randrange, choice
import uuid
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from final import RecordCollection, read_fasta_file

def generate_sequence(sequence_length):
    return ''.join(choice("ACGT") for _ in range(sequence_length))    

def create_record(index, sequence_length = None, id = None):
    sequence = generate_sequence(randrange(500, 5000) if not sequence_length else sequence_length)
    record = SeqRecord(Seq(sequence), id = str(uuid.uuid4()) if not id else id)
    record.record_number = index
    return record

def generate_records(num_records):
    return [create_record(i) for i in range(num_records)]

class SequenceTestCase(unittest.TestCase):
    def test_fasta_read(self):
        records = read_fasta_file("dna.example.fasta")
        self.assertEqual(len(records), 25)
        for i, record in enumerate(records):
            self.assertEqual(i, record.record_number)
            self.assertIsInstance(record, SeqRecord)

    def test_get_record_by_id(self):
        records = []
        ids = []
        fake_ids = []
        num_records = 50

        for i in range(num_records):
            id = str(uuid.uuid4())
            records.append(create_record(i, sequence_length = None, id = id))
            ids.append(id)
            fake_ids.append(str(uuid.uuid4()))

        records = RecordCollection(records)

        for i in range(num_records):
            self.assertIsInstance(records.get_record_by_id(ids[i]), SeqRecord)
            self.assertIsNone(records.get_record_by_id(fake_ids[i]))

    
    def test_longest_and_shortest_sequence_lengths(self):
        records = generate_records(50)

        i = randrange(len(records) + 1)
        max_length_record = create_record(i, sequence_length = 5001)
        records.insert(i, max_length_record)

        i = randrange(len(records) + 1)
        min_length_record = create_record(i, sequence_length = 499)
        records.insert(i, min_length_record)

        records = RecordCollection(records)
        longest_records = records.get_longest_sequence_records()
        shortest_records = records.get_shortest_sequence_records()
        
        self.assertGreater(len(records), 0)
        self.assertEqual(max_length_record.id, longest_records[0].id)
        self.assertEqual(min_length_record.id, shortest_records[0].id)    

    def test_orf(self):
        records = read_fasta_file("dna2.fasta")
        records = RecordCollection(records)

if __name__ == '__main__':
    unittest.main()