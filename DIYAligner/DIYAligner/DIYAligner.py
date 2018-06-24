class Aligner():

    def __init__(self):
        self.complements = {'A' : 'T', 
                            'C' : 'G', 
                            'G' : 'C', 
                            'T' : 'A'}

    def parse_fasta(self, filename):
        reference = ""
        with open(filename, "r") as file:
            file.readline()
            reference = file.readline()
        return reference

    def parse_fastq(self, filename):
        reads = []
        with open(filename, "r") as file:
            while True:
                 first_line = f.readline()
                 if len(first_line) == 0:
                     break
                 name = first_line[1:].rstrip()
                 seq = f.readline().rstrip()
                 f.readline()
                 qual = f.readline().rstrip()
                 reads.append((name, seq, qual))
        return reads

    def complement_reads(self, reads):
        reads_comp = []
        for read in reads:
            comp_read = ""
            for letter in read:
                comp_read += self.complements[letter]
            reads_comp.append((read[0], comp_read, read[2]))
        return reads_comp
