import sys
import numpy
import statistics

class Aligner():

    def __init__(self):
        self.complements = {'A' : 'T', 
                            'C' : 'G', 
                            'G' : 'C', 
                            'T' : 'A'}
        self.field_size = -1

    def parse_fasta(self, filename):
        self.reference = ""
        with open(filename, "r") as file:
            file.readline()
            self.reference = file.readline()

    def parse_fastq(self, filename):
        self.reads = []
        with open(filename, "r") as f:
            while True:
                 first_line = f.readline()
                 if len(first_line) == 0:
                     break
                 name = first_line[1:].rstrip()
                 seq = f.readline().rstrip()
                 f.readline()
                 qual = f.readline().rstrip()
                 self.reads.append((name, seq, qual))

    def complement_reads(self):
        self.reads_comp = []
        for read in self.reads:
            comp_read = ""
            for letter in read[1]:
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
        seeds = set()
        pos = 0
        while True:
            if len(read) - pos < length:
                break
            seeds.add(read[pos:pos+length])
            pos += interval
        return seeds  

    def seed(self, read, length, interval):
        if self.field_size == -1 or self.field_size != length:
            self.create_index(length)
        seeds = self.create_seeds(read, length, interval)
        seed_hits = {}
        for seed in seeds:
            seed_hits[seed] = self.index[seed]
        return seed_hits

    def choose_extend_place(self, seed_hits, read_length):
        all_hits =[]
        for seed in seed_hits:
            all_hits.extend(seed_hits[seed])
        hits_in_area = {}
        for hit in all_hits:
            hits_in_area[hit] = [hit]
            for i in all_hits:
                if abs(hit - i) < read_length and hit != i:
                    hits_in_area[hit].append(i)
        max_neighbours = []
        for hit in hits_in_area:
            if len(hits_in_area[hit]) >  len(max_neighbours):
                max_neighbours = hits_in_area[hit]
        extend_location = statistics.median(max_neighbours)
        return extend_location
    
    def scoring_matrix(self, a, b):
        if a==b: return 1 
        if a=='_' or b=='_': return -7 
        return -4 

    def local_alignment(self, x, y, s):
        D = numpy.zeros((len(x) + 1,len(y) + 1), dtype=int)

        for i in range(1,len(x) + 1):
            for j in range(1,len(y) + 1):
                D[i,j] = max(0, 
                       D[i - 1,j] + s(x[i - 1], '_'),
                       D[i,j - 1] + s('_',   y[j - 1]), 
                       D[i - 1,j - 1] + s(x[i - 1],y[j - 1]))
    
        local_max = D.max()
        return D, local_max

    def tracebak(self, x, y, V, s):
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

    def align_read(self, read, length = 20, interval = 10):
        seed_hits = self.seed(read, length, interval)
        extend_location = self.choose_extend_place(seed_hits, len(read))
        ref = self.reference[int(extend_location - len(read) - 1):int(extend_location + len(read) + 1)]
        D, max = self.local_alignment(ref, read, self.scoring_matrix)
        aligment, transcript = self.tracebak(ref, read, D, self.scoring_matrix)
        print(aligment)
        print(transcript)

al = Aligner()
al.parse_fasta("InputFiles/MT.fa.txt")
#print(len(al.reference)) 16570
al.parse_fastq("InputFiles/synth.fq.txt")
#for read in al.reads:
#    print('Name = {}, Sequence = {}, Quality = {} \n'.format(read[0], read[1], read[2]))
#print(len(al.reads)) 2488
al.complement_reads()
#print(al.reads[0][1])
#print(al.reads_comp[0][1][::-1])
#al.create_index(20)
#for i in range(10000):
#    print(al.index[al.reference[i:20 + i]])
#print(al.reads[0][1])
#seeds = al.create_seeds(al.reads[0][1], 20, 10))
#seed_hits = al.seed(al.reads[0][1], 20, 10)
#for hit in seed_hits:
#    print(seed_hits[hit])
#extend_place = al.choose_extend_place(seed_hits, len(al.reads[0][1]))
#print(extend_place)
al.align_read(al.reads[0][1])




