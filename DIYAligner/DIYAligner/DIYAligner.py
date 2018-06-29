import sys
import numpy as np

class Aligner():

    def __init__(self):
        self.complements = {'A' : 'T', 
                            'C' : 'G', 
                            'G' : 'C', 
                            'T' : 'A'}
        self.field_size = -1

    def parse_fasta(self, filename, field_length):
        self.reference = ""
        with open(filename, "r") as file:
            file.readline()
            self.reference = file.readline()

    def parse_fastq(self, filename):
        self.reads = []
        with open(filename, "r") as file:
            while True:
                 first_line = f.readline()
                 if len(first_line) == 0:
                     break
                 name = first_line[1:].rstrip()
                 seq = f.readline().rstrip()
                 f.readline()
                 qual = f.readline().rstrip()
                 self.reads.append((name, seq, qual))

    def complement_reads(self, reads):
        self.reads_comp = []
        for read in reads:
            comp_read = ""
            for letter in read:
                comp_read = self.complements[letter] + comp_read
            self.reads_comp.append((read[0], comp_read, read[2]))

    def create_index(self, ln):
        self.index = {}
        self.field_size = ln
        for i in range(len(self.reference) - ln + 1):
            substr = self.reference[i:i+ln]
            if substr in self.index:
                self.index[substr].append(i)
            else:
                self.index[substr] = [i]

    def create_seeds(self, read, length, interval):
        seeds = []
        pos = 0
        while True:
            if len(read) - pos < length:
                break
            seeds.append(read[i:i+length])
            pos += interval
        return seeds  

    def seed(self, read, lenght = 20, interval = 10):
        if self.field_size == -1 or self.field_size != length:
            create_index(length)
        seeds = create_seeds(read, lenght, interval)
        seed_hits = []
        for seed in seeds:
            seed_hits.append(self.index[seed])
        

