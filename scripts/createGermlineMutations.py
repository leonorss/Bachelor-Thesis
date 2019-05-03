import vcf
import random
from Bio import SeqIO

# read in parameters from config File
n_g_m = snakemake.params[0]
f_h_m = snakemake.params[1]

# create a seed for reproducibilaty
configSeed = snakemake.params[2]
random.seed(configSeed)

# number if Homo- and Heterozygous Mutations according to the config File
numberOfHeterozygousGermlineMutations = int(n_g_m * f_h_m)
numberOfHomozygousGermlineMutations = n_g_m - numberOfHeterozygousGermlineMutations

# read in reference Genome to check that we don't insert Mutations, that actually don't change the Genome
referenceGenome = SeqIO.read(snakemake.input[0], "fasta")
referenceGenomeLength = len(referenceGenome)

# create two vcf Files for the 2 allels with the given header
vcfReader_MetaData = vcf.Reader(filename=snakemake.input[1])
vcfWriter1 = vcf.Writer(open(snakemake.output[0], 'w'), vcfReader_MetaData)
vcfWriter2 = vcf.Writer(open(snakemake.output[1], 'w'), vcfReader_MetaData)

# creating specified number of random homozygous mutations: 1=A, 2=C, 3=G, 4=T
for homozygousMutation in range(0, numberOfHomozygousGermlineMutations):
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
                if homozygousMutation != 0:
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                # inserting the mutations into the two vcf files
                if alreadyMutated == 0:
                    inserted = 1
                    record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                                REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                                ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                                FORMAT=".", sample_indexes=[], samples=None)
                    vcfWriter1.write_record(record)
                    vcfWriter2.write_record(record)
                    vcfWriter1.flush()
                    vcfWriter2.flush()

        # case nucleotid="C"
        elif newMutation==2:
            mutatedNucleotid = "C"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("c" != referenceGenome[insertPlace]):
                if homozygousMutation != 0:
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break

                # inserting the mutations into the two vcf files
                if alreadyMutated == 0:
                    inserted = 1
                    record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                                REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                                ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                                FORMAT=".", sample_indexes=[], samples=None)
                    vcfWriter1.write_record(record)
                    vcfWriter2.write_record(record)
                    vcfWriter1.flush()
                    vcfWriter2.flush()


        # case nucleotid="G"
        elif newMutation==3:
            mutatedNucleotid = "G"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("g" != referenceGenome[insertPlace]):
                if homozygousMutation != 0:
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break

                # inserting the mutations into the two vcf files
                if alreadyMutated == 0:
                    inserted = 1
                    record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                                REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                                ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                                FORMAT=".", sample_indexes=[], samples=None)
                    vcfWriter1.write_record(record)
                    vcfWriter2.write_record(record)
                    vcfWriter1.flush()
                    vcfWriter2.flush()

        # case nucleotid="T"
        else:
            mutatedNucleotid = "T"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("t" != referenceGenome[insertPlace]):
                if homozygousMutation != 0:
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break

                # inserting the mutations into the two vcf files
                if alreadyMutated == 0:
                    inserted = 1
                    record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                                REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                                ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                                FORMAT=".", sample_indexes=[], samples=None)
                    vcfWriter1.write_record(record)
                    vcfWriter2.write_record(record)
                    vcfWriter1.flush()
                    vcfWriter2.flush()


# creating specified number of random heterozygous mutations: 1=A, 2=C, 3=G, 4=T
for heterozygousMutation in range(0, numberOfHeterozygousGermlineMutations):
    inserted = 0

    # choose either Allel randomly
    allel = random.randint(1, 2)

    # checking that we didn't create a mutation that din't actually change the reference genome
    while(inserted == 0):

        # choosing random Nucleotid and random place to insert
        newMutation = random.randrange(1, 4)
        insertPlace = random.randrange(0, (referenceGenomeLength-1))
        alreadyMutated = 0

        # case nucleotid="A"
        if newMutation==1:
            mutatedNucleotid = "A"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("a" != referenceGenome[insertPlace]):
                record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                            REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                            ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                            FORMAT=".", sample_indexes=[], samples=None)

                # insert record into randomly choosen allel
                if (allel == 1):
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter1.write_record(record)
                        vcfWriter1.flush()
                else:
                    vcfReader2 = vcf.Reader(filename=snakemake.output[1])
                    for existingRecord in vcfReader2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter2.write_record(record)
                        vcfWriter2.flush()

        # case nucleotid="C"
        elif newMutation==2:
            mutatedNucleotid = "C"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("c" != referenceGenome[insertPlace]):
                record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                            REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                            ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                            FORMAT=".", sample_indexes=[], samples=None)

                # insert record into randomly choosen allel
                if (allel == 1):
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter1.write_record(record)
                        vcfWriter1.flush()
                else:
                    vcfReader2 = vcf.Reader(filename=snakemake.output[1])
                    for existingRecord in vcfReader2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter2.write_record(record)
                        vcfWriter2.flush()

        # case nucleotid="G"
        elif newMutation==3:
            mutatedNucleotid = "G"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("g" != referenceGenome[insertPlace]):
                record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                            REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                            ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                            FORMAT=".", sample_indexes=[], samples=None)

                # insert record into randomly choosen allel
                if (allel == 1):
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter1.write_record(record)
                        vcfWriter1.flush()
                else:
                    vcfReader2 = vcf.Reader(filename=snakemake.output[1])
                    for existingRecord in vcfReader2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter2.write_record(record)
                        vcfWriter2.flush()

        # case nucleotid="T"
        else:
            mutatedNucleotid = "T"

            # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
            if (mutatedNucleotid != referenceGenome[insertPlace]) and ("t" != referenceGenome[insertPlace]):
                record = vcf.model._Record(CHROM=22, POS=(insertPlace+1), ID='.',
                            REF=vcf.model._Substitution(referenceGenome[insertPlace]),
                            ALT=[vcf.model._Substitution(mutatedNucleotid)], QUAL='.', FILTER='PASS', INFO={},
                            FORMAT=".", sample_indexes=[], samples=None)

                # insert record into randomly choosen allel
                if (allel == 1):
                    vcfReader1 = vcf.Reader(filename=snakemake.output[0])
                    for existingRecord in vcfReader1:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter1.write_record(record)
                        vcfWriter1.flush()
                else:
                    vcfReader2 = vcf.Reader(filename=snakemake.output[1])
                    for existingRecord in vcfReader2:
                        if insertPlace == existingRecord.POS:
                            alreadyMutated = 1
                            break
                    if alreadyMutated == 0:
                        inserted = 1
                        vcfWriter2.write_record(record)
                        vcfWriter2.flush()

# close writers
vcfWriter1.close()
vcfWriter2.close()







# f1= open(snakemake.output[0],"w+")
# f2= open(snakemake.output[1],"w+")
# f1.write("%i" %(n_g_m))
# f2.write("%i" %(n_g_m))
# f1.close()
# f2.close()
