import sys
import numpy as np
import statistics

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

    def seed(self, read, lenght, interval):
        if self.field_size == -1 or self.field_size != length:
            create_index(length)
        seeds = create_seeds(read, lenght, interval)
        seed_hits = {}
        for seed in seeds:
            seed_hits[seed] = self.index[seed]
        return seed_hits


    def choose_extend_place(self, seed_hits, read_length):
        all_hits =[]
        for seed in seed_hits:
            all_hits.append(seed_hits[seed])
        hits_in_area = {}
        for hit in seed_hits:
            hits_in_area[hit] = [hit]
            for i in seed_hits:
                if abs(hit - i) < read_length and hit != i:
                    hits_in_area[hit].append[i]
        max_neighbours = -1
        for hit in hits_in_area:
            if len(hits_in_area[hit]) >  max_neighbours:
                max_neighbours = hits_in_area[hit]
        max_neighbours = max_neighbours.sort()
        extend_location = statistics.median(max_neighbours)
        return extend_location
       

