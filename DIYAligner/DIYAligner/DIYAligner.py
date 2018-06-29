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
    
    def scoring_matrix(a,b):
        if a==b: return 1 
        if a=='_' or b=='_': return -7 
        return -4 

    def local_alignment(x, y, s):
        D = numpy.zeros((len(x) + 1,len(y) + 1), dtype=int)

        for i in range(1,len(x) + 1):
            for j in range(1,len(y) + 1):
                D[i,j] = max(0, 
                       D[i - 1,j] + s(x[i - 1], '_'),
                       D[i,j - 1] + s('_',   y[j - 1]), 
                       D[i - 1,j - 1] + s(x[i - 1],y[j - 1]))
    
        local_max = D.max()
        return D, local_max

    def tracebak(x,y,V, s):
        maxValue=numpy.where(V==V.max())
        i=maxValue[0][0]
        j=maxValue[1][0]
        ax,ay,am, tr = '','','',''
        while (i>0 or j>0) and V[i,j]!=0:
            if i>0 and j>0:
                sigma = 1 if x[i-1]==y[j-1] else 0
                d=V[i-1,j-1]+s(x[i-1],y[j-1]) # diagonal move
            if i>0: v=V[i-1,j]+s(x[i-1],'_')  # vertical move
            if j>0: h=V[i,j-1]+s('_',y[j-1])  # horizontal move
        
            # diagonal is the best
            if d>=v and d>=h:
                ax += x[i-1]
                ay += y[j-1]
                if sigma==1:
                    tr+='M'
                    am+='|'
                else:
                    tr+='R'
                    am+=' '
                i-=1
                j-=1
            # vertical is the best
            elif v>=h:
                ax+=x[i-1]
                ay+='_'
                tr+='I'
                am+=' '
                i-=1
            # horizontal is the best
            else:
                ay+=y[j-1]
                ax+='_'
                tr+='D'
                am+=' '
                j-=1
        alignment= '\n'.join([ax[::-1], am[::-1], ay[::-1]])
        return alignment, tr[::-1]


