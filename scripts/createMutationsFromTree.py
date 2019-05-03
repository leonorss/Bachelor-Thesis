import vcf
import random
from Bio import SeqIO
import numpy as np

# create a seed for reproducibilaty
configSeed = snakemake.params[2]
random.seed(configSeed)

# number of mutations: if one mutation at every internal node except the root n-2
numberOfMutations = snakemake.params[0]

# read in reference Genome to check that we don't insert Mutations, that actually don't change the Genome
referenceGenome = SeqIO.read(snakemake.input[3], "fasta")
referenceGenomeLength = len(referenceGenome)

# read in tree and number of leaves
tree = np.loadtxt(fname = snakemake.input[0])
leaves = snakemake.params[1]

# read in metadata for the vcf files
vcfReader_MetaData = vcf.Reader(filename=snakemake.input[4])
vcfReader_Allel1 = vcf.Reader(filename=snakemake.input[1])
vcfReader_Allel2 = vcf.Reader(filename=snakemake.input[2])

# save created mutations and positions to insert into vcf files later
randomPostions = np.zeros(numberOfMutations)
mutatedNucleotids = np.chararray(numberOfMutations)

for mutation in range(0, numberOfMutations):
    inserted = 0

    # checking that we didn't create a mutation that din't actually change the reference genome
    while(inserted == 0):

        # choosing random Nucleotid and random place to insert
        newMutation = random.randrange(1, 5)
        insertPlace = random.randrange(0, referenceGenomeLength)
        alreadyMutated = 0

        # case nucleotid="A"
        if newMutation==1:
            mutatedNucleotid = "A"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("a" != referenceGenome[insertPlace]):

                # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                vcfReader_Allel1 = vcf.Reader(filename=snakemake.input[1])

                for existingRecord in vcfReader_Allel1:

                    if insertPlace == existingRecord.POS:
                        alreadyMutated = 1
                        break

                if alreadyMutated != 1:

                    # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                    vcfReader_Allel2 = vcf.Reader(filename=snakemake.input[2])

                    for existingRecord in vcfReader_Allel2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break

                # inserting the mutations and their insert positions into the arrays randomPostions and mutatedNucleotids
                if alreadyMutated == 0:
                    inserted = 1
                    randomPostions[mutation] = insertPlace
                    mutatedNucleotids[mutation] = mutatedNucleotid

        # case nucleotid="C"
        elif newMutation==2:
            mutatedNucleotid = "C"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("c" != referenceGenome[insertPlace]):

                # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                vcfReader_Allel1 = vcf.Reader(filename=snakemake.input[1])

                for existingRecord in vcfReader_Allel1:
                    if insertPlace == existingRecord.POS:
                        alreadyMutated = 1
                        break

                if alreadyMutated != 1:

                    # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                    vcfReader_Allel2 = vcf.Reader(filename=snakemake.input[2])

                    for existingRecord in vcfReader_Allel2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break

                # inserting the mutations and their insert positions into the arrays randomPostions and mutatedNucleotids
                if alreadyMutated == 0:
                    inserted = 1
                    randomPostions[mutation] = insertPlace
                    mutatedNucleotids[mutation] = mutatedNucleotid

        # case nucleotid="G"
        elif newMutation==3:
            mutatedNucleotid = "G"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("g" != referenceGenome[insertPlace]):

                # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                vcfReader_Allel1 = vcf.Reader(filename=snakemake.input[1])

                for existingRecord in vcfReader_Allel1:
                    if insertPlace == existingRecord.POS:
                        alreadyMutated = 1
                        break

                if alreadyMutated != 1:

                    # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                    vcfReader_Allel2 = vcf.Reader(filename=snakemake.input[2])

                    for existingRecord in vcfReader_Allel2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break

                # inserting the mutations and their insert positions into the arrays randomPostions and mutatedNucleotids
                if alreadyMutated == 0:
                    inserted = 1
                    randomPostions[mutation] = insertPlace
                    mutatedNucleotids[mutation] = mutatedNucleotid

        # case nucleotid="T"
        else:
            mutatedNucleotid = "T"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("t" != referenceGenome[insertPlace]):

                # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                vcfReader_Allel1 = vcf.Reader(filename=snakemake.input[1])

                for existingRecord in vcfReader_Allel1:
                    if insertPlace == existingRecord.POS:
                        alreadyMutated = 1
                        break

                if alreadyMutated != 1:

                    # check that we don't mutate the same nucleotid as in we already did in the germline mutations
                    vcfReader_Allel2 = vcf.Reader(filename=snakemake.input[2])

                    for existingRecord in vcfReader_Allel2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break

                # inserting the mutations and their insert positions into the arrays randomPostions and mutatedNucleotids
                if alreadyMutated == 0:
                    inserted = 1
                    randomPostions[mutation] = insertPlace
                    mutatedNucleotids[mutation] = mutatedNucleotid

# for every leave saving its mutations: row=mutations, columns=leave
mutationsOfAllLeaves = np.zeros([numberOfMutations, (2*leaves)])

# insert mutations of a node into all it's children
def insertMutationIntoChildren(parent, mutation):

    if (tree[2][parent] == -1):
        mutationsOfAllLeaves[mutation][int(tree[0][parent])] = 1
        mutationsOfAllLeaves[mutation][int(tree[1][parent])] = 1

        insertMutationIntoChildren(int(tree[0][parent]), mutation)
        insertMutationIntoChildren(int(tree[1][parent]), mutation)

# randomly distribute mutations in the tree
for mutation in range(0, numberOfMutations):

    inserted = 0
    while (inserted == 0):

        insertionNode = random.randrange(1, int(tree.size/3))

        if (tree[2][insertionNode] == (-1)):
            mutationsOfAllLeaves[mutation][insertionNode] = 1
            insertMutationIntoChildren(insertionNode, mutation)
            inserted = 1

# create output files
for leave in range(0, leaves):

    vcf_writer1 = vcf.Writer(open(snakemake.output[leave], "w"), vcfReader_MetaData)
    vcf_writer2 = vcf.Writer(open(snakemake.output[leaves + leave], "w"), vcfReader_MetaData)

    vcfReader_Allel1 = vcf.Reader(filename=snakemake.input[1])
    for record in vcfReader_Allel1:
        vcf_writer1.write_record(record)
        vcf_writer1.flush()

    vcfReader_Allel2 = vcf.Reader(filename=snakemake.input[2])
    for record in vcfReader_Allel2:
        vcf_writer2.write_record(record)
        vcf_writer2.flush()

    for mutation in range(0, numberOfMutations):

        if (mutationsOfAllLeaves[mutation][leave] == 1):
            position = int(randomPostions[mutation])
            mutation = mutatedNucleotids[mutation]
            record = vcf.model._Record(CHROM=22, POS=(position+1), ID='.',
                        REF=vcf.model._Substitution(referenceGenome[position]),
                        ALT=[vcf.model._Substitution(mutation)], QUAL='.', FILTER='PASS', INFO={},
                        FORMAT=".", sample_indexes=[], samples=None)
            vcf_writer1.write_record(record)
            vcf_writer2.write_record(record)
            vcf_writer1.flush()
            vcf_writer2.flush()

    vcf_writer1.close()
    vcf_writer2.close()
