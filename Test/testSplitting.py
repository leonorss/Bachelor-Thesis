import numpy as np
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

random.seed(234)

genome = SeqIO.read("ReferenceGenome.fa", "fasta")
genomeLength = len(genome)

mean = 5
standardDeviation = 0.9

reads = [0]

i = 0
position = 0
while (position < (genomeLength - 2*mean)):
    length = np.random.normal(mean, standardDeviation)
    length = int(np.round(length))

    position = length + reads[i] - 1
    reads.append(position)
    i = i + 1

length = ((genomeLength - 1) - reads[i]) / 2
length = int(np.round(length))
position = length + reads[i]
reads.append(position)

reads.append(genomeLength)

print(reads)

i=0
for read in range(0, (len(reads)-1)):
    seq_Record = SeqRecord((genome[(reads[read]):(reads[read+1])].seq), id=(genome.id))
    name = "Record_" + str(i) + ".fa"
    SeqIO.write(seq_Record, name, "fasta")
    i = i + 1
