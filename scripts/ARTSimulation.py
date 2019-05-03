import numpy as np
import random
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord

# read in seed from config file for reproducibilaty
configSeed = snakemake.params[0] + snakemake.wildcards.sample
random.seed(configSeed)

# read in reference genome
referenceGenome = SeqIO.read(snakemake.input[1], "fasta")

# read in mean and variance of segment length from config file
meanBinSizeLength = snakemake.params[1]
standardDeviationOfBinSize = snakemake.params[2]

# read in genome sequence
genome = SeqIO.read(snakemake.input, "fasta")
genomeLength = len(genome)

# read in mean and standart deviation of coverage from config file
meanCoverage = snakemake.params[3]
standardDeviationCoverage = snakemake.params[4]

# read in read length from config file
readLength = snakemake.params[5]

# read in mean and standart deviation of fragment size from config file
meanFragmentSize = snakemake.params[6]
standardDeviationOfFragmentSize = snakemake.params[7]

# read in the probability of an MDA amplification error
MDAamplificationErrorProbability = snakemake.params[8]

# save end positions of read lengths in vector
bins = [0]

# randomly generate bin lengths with given parameters and save the end positions with respect to the whole genome into the bins vector
i = 0
position = 0
while (position < (genomeLength - 2*meanBinSize)):
    length = np.random.normal(meanBinSize, standardDeviationOfBinSize)
    length = np.round(length)

    position = length + bins[i]
    bins.append(position)
    i = i + 1

# to make sure we don't have a very small bin at the end because of standart deviation, we calculate the last 2 bin sizes,
# by stopping the above loop at (genomeLength - 2*meanBinSize) and dividing the rest genome by 2
length = (genomeLength - reads[i]) / 2
length = np.round(length)
position = length + bins[i]
bins.append(position)

bins.append(genomeLength)

# now we create the according fasta files
for bin in range(0, (len(bins)-1)):
    seq_Record = SeqRecord((genome[(bins[bin]):(bins[bin+1])].seq), id=(genome.id))
    nametemplate = "Bin_Tree" + str(snakemake.wildcards.treename) + "_Sample" + str(snakemake.wildcards.sample) + "_Bin" + str(i) + "_"
    name = nametemplate + "0.fa"
    SeqIO.write(seq_Record, name, "fasta")

    readCoverage = np.random.normal(meanCoverage, standardDeviationCoverage)
    readCoverage = np.round(readCoverage)

    binLength = bins[bin+1] - bins[bin]
    for amplification in range(1, readCoverage):
        choosenReadNumber = random.randrange(1, (amplification+1))
        amplificationError = random.choices(population = [0, 1], weights = [(1-MDAamplificationErrorProbability), MDAamplificationErrorProbability], k = 1)
        choosenReadName = nametemplate + str(choosenReadNumber) + ".fa"
        newName = nametemplate + str(amplification) + ".fa"

        choosenRead = SeqIO.read(choosenReadName, "fasta")
        seq_Record = SeqRecord(choosenRead.seq, id=(choosenRead.id))
        SeqIO.write(seq_Record, newName, "fasta")

        if (amplificationError == 1):
            inserted = 0

            # checking that we didn't create a mutation that din't actually change the reference genome
            while(inserted == 0):
                # choosing random Nucleotid and random place to insert
                newMutation = random.randrange(1, 5)
                insertPlace = random.randrange(0, binLength)
                alreadyMutated = 0

                # case nucleotid="A"
                if newMutation==1:
                    mutatedNucleotid = "A"

                    # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
                    if (mutatedNucleotid != referenceGenome[insertPlace]) and ("a" != referenceGenome[insertPlace]):
                        vcfReader = vcf.Reader(filename=choosenReadName)
                        for existingRecord in vcfReader:
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
                        vcfWriter.write_record(record)
                        vcfWriter.flush()


                # case nucleotid="C"
                elif newMutation==2:
                    mutatedNucleotid = "C"

                    # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
                    if (mutatedNucleotid != referenceGenome[insertPlace]) and ("c" != referenceGenome[insertPlace]):

                        vcfReader = vcf.Reader(filename=choosenReadName)
                        for existingRecord in vcfReader:
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
                        vcfWriter.write_record(record)
                        vcfWriter.flush()



                # case nucleotid="G"
                elif newMutation==3:
                    mutatedNucleotid = "G"

                    # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
                    if (mutatedNucleotid != referenceGenome[insertPlace]) and ("g" != referenceGenome[insertPlace]):
                        vcfReader = vcf.Reader(filename=choosenReadName)
                        for existingRecord in vcfReader:
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
                        vcfWriter.write_record(record)
                        vcfWriter.flush()


                # case nucleotid="T"
                else:
                    mutatedNucleotid = "T"

                    # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
                    if (mutatedNucleotid != referenceGenome[insertPlace]) and ("t" != referenceGenome[insertPlace]):
                        vcfReader = vcf.Reader(filename=choosenReadName)
                        for existingRecord in vcfReader:
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
                        vcfWriter.write_record(record)
                        vcfWriter.flush()

            # close writers
            vcfWriter.close()

        outputfile = "ARTRecord_Tree" + str(snakemake.wildcards.treename) + "_Sample" + str(snakemake.wildcards.sample) + "_Bin" + str(i) + "_" + str(amplification)
        shellCommand = "art_illumina -na -i " + newName + " -p -l " + str(readLength) + " -ss HS25 -f 3 -m " + str(meanFragmentSize) + " -s " + str(standardDeviationOfFragmentSize) + " -o " + outputfile
        os.system(shellCommand)
