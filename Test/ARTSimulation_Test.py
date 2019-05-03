import numpy as np
import random
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord

# genome = SeqIO.read("ReferenceGenome.fa", "fasta")
#
# seq_Record = SeqRecord((genome[3:5].seq), id=(genome.id))
# name = "Record_" + str(1) + ".fa"
# SeqIO.write(seq_Record, name, "fasta")
#
# readCoverage = 5
# outputfile = "ARTRecord_" + str(1)
# shellCommand = "art_illumina -sam -na -i " + name + " -l " + str(2) + " -ss HS25 -f " + str(readCoverage) +  " -o " + outputfile
# os.system(shellCommand)

prob = 0.1
# x = np.random.choice(2, size=20, p=[(1-prob), prob])
# print(x)

x = random.choices(population=[0, 1], weights=[(1-prob), prob], k=10)
print(x)
