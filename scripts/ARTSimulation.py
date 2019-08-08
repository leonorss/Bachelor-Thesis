import numpy as np
import random
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import vcf
import time

# create necessary folders
os.system("mkdir results")
os.system("mkdir Data/MDASimulation")

# create name for the vcf file
vcfName = "Data/MDASimulation/MDAerrors_" + snakemake.wildcards.treename + "_Allel" + snakemake.wildcards.sample + ".vcf"

# save created mutations and positions to insert into vcf file later
randomPostions = []
mutatedNucleotids = []

# read in reference Genome to get the reference nucleotid
referenceGenome = SeqIO.read(snakemake.params[14], "fasta")

# read in seed from config file for reproducibilaty
configSeed = snakemake.params[0] + int(snakemake.wildcards.sample)
random.seed(configSeed)

# read in mean and variance of segment length from config file
meanBinSizeLength = snakemake.params[1]
standartDeviationOfBinSize = snakemake.params[2]
varianceOfBinSize = standartDeviationOfBinSize * standartDeviationOfBinSize

# read in genome sequence
genome = SeqIO.read(snakemake.input[0], "fasta")
genomeLength = len(genome)

# read in mean and standart deviation of coverage from config file
meanCoverage = snakemake.params[3]
standardDeviationOfCoverage = snakemake.params[4]
varianceOfCoverage = standardDeviationOfCoverage * standardDeviationOfCoverage
zeroInflatationProb = snakemake.params[12]

# read in read length from config file
readLength = snakemake.params[5]

# read in mean and standart deviation of fragment size from config file
meanFragmentSize = snakemake.params[6]
standartDeviationOfFragmentSize = snakemake.params[7]
varianceOfFragmentSize= standartDeviationOfFragmentSize * standartDeviationOfFragmentSize

# read in the probability of an MDA amplification error
MDAamplificationErrorProbability = snakemake.params[8]

# save the name of the choosen Chromosom
chromosomReader = vcf.Reader(filename=snakemake.params[9])
chromosom = int(next(iter(chromosomReader.contigs)))

# save the ART Illumina coverage
ARTcoverage = snakemake.params[10]

# save the probability of a read being dropped during sequencing
ARTDropReadProbability = snakemake.params[11]

# save end positions of read lengths in vector
bins = [0]

# randomly generate bin lengths with given parameters and save the end positions with respect to the whole genome into the bins vector
i = 0
position = 0
while (position < (genomeLength - 2*meanBinSizeLength)):

    # generating a random bin size with a negative binomial distribution and parameters p and r
    rBinSize = (meanBinSizeLength*meanBinSizeLength) / (varianceOfBinSize - meanBinSizeLength)
    pBinSize = (varianceOfBinSize - meanBinSizeLength) / varianceOfBinSize
    length = np.random.negative_binomial(rBinSize, (1 - pBinSize))

    position = length + bins[i]
    bins.append(position)
    i = i + 1

# to make sure we don't have a very small bin at the end because of standart deviation, we calculate the last 2 bin sizes,
# by stopping the above loop at (genomeLength - 2*meanBinSize) and dividing the rest genome by 2
length = (genomeLength - bins[i]) / 2
length = int(np.round(length))
position = length + bins[i]
bins.append(position)

bins.append(genomeLength)

# creating two result vectors to save all the FASTQ outputs and create two paired ends read files later
ArtIlluminaRecords1 = []
ArtIlluminaRecords2 = []

# now we create the according fasta files
for bin in range(0, (len(bins)-1)):

    # generating a random read coverage with a negative binomial distribution and parameters p and r
    rReadCoverage = (meanCoverage*meanCoverage) / (varianceOfCoverage - meanCoverage)
    pReadCoverage = (varianceOfCoverage - meanCoverage) / varianceOfCoverage

    whichDistribution = random.choices(population = [0, 1], weights = [zeroInflatationProb, (1-zeroInflatationProb)], k = 1)

    if whichDistribution[0] == 0:
        readCoverage = 0
    else:
        readCoverage = np.random.negative_binomial(rReadCoverage, (1 - pReadCoverage))

    if (readCoverage != 0):

        # sequence of the bin according to the length specified in the bins vector
        seqRecord = SeqRecord((genome[(bins[bin]):(bins[bin+1])].seq), id=(str(bin) + "." + str(0)))
        nametemplate = "Bin_Tree" + snakemake.wildcards.treename + "_Sample" + snakemake.wildcards.sample + "_Bin" + str(bin) + ".fa"

        # creating a first record with the sequence above and save it in the records vector
        records = [seqRecord]

        binLength = bins[bin+1] - bins[bin]

        # to simulate the MDA simulation we amplificate the sequence according to the choosen coverage and introduce amlification errors
        for amplification in range(1, readCoverage):
            # randomly choose a sequence generated by the amplification to be copied
            choosenReadNumber = random.randrange(0, (amplification))
            # randomly choose if a MDA amplification error gets introduced or not with a probability= MDAamplificationErrorProbability
            amplificationError = random.choices(population = [0, 1], weights = [(1-MDAamplificationErrorProbability), MDAamplificationErrorProbability], k = 1)

            # copy choosen sequence
            choosenReadSeq = (records[choosenReadNumber]).seq
            idName = str(bin) + "." + str(amplification)

            # if there's a MDA amlification error ocurred, we create it and introduce it into the copy of the choosen sequence
            if (amplificationError[0] == 1):
                inserted = 0

                # checking that we didn't create a mutation that didn't actually change the reference genome
                while(inserted == 0):
                    # choosing random Nucleotid and random place to insert
                    newMutation = random.randrange(1, 5)
                    insertPlace = random.randrange(0, binLength)

                    # case nucleotid="A"
                    if newMutation==1:
                        mutatedNucleotid = "A"

                        # making sure the nucleotid is not the same as in the given sequence or can not be sequenced
                        if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("a" != choosenReadSeq[insertPlace]) and ("N" != choosenReadSeq[insertPlace]):

                            # inserting the mutations into the copied sequence
                            mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                            mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                            newSeq = mutableChoosenReadSeq.toseq()

                            newRecord = SeqRecord(newSeq, id=idName)

                            randomPostions.append(bins[bin]+insertPlace+1)
                            mutatedNucleotids.append(mutatedNucleotid)

                            inserted = 1




                    # case nucleotid="C"
                    elif newMutation==2:
                        mutatedNucleotid = "C"

                        # making sure the nucleotid is not the same as in the given sequence or can not be sequenced
                        if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("c" != choosenReadSeq[insertPlace]) and ("N" != choosenReadSeq[insertPlace]):

                            # inserting the mutations into the copied sequence
                            mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                            mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                            newSeq = mutableChoosenReadSeq.toseq()

                            newRecord = SeqRecord(newSeq, id=idName)

                            randomPostions.append(bins[bin]+insertPlace+1)
                            mutatedNucleotids.append(mutatedNucleotid)

                            inserted = 1



                    # case nucleotid="G"
                    elif newMutation==3:
                        mutatedNucleotid = "G"

                        # making sure the nucleotid is not the same as in the given sequence or can not be sequenced
                        if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("g" != choosenReadSeq[insertPlace]) and ("N" != choosenReadSeq[insertPlace]):

                            # inserting the mutations into the copied sequence
                            mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                            mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                            newSeq = mutableChoosenReadSeq.toseq()

                            newRecord = SeqRecord(newSeq, id=idName)

                            randomPostions.append(bins[bin]+insertPlace+1)
                            mutatedNucleotids.append(mutatedNucleotid)

                            inserted = 1


                    # case nucleotid="T"
                    else:
                        mutatedNucleotid = "T"

                        # making sure the nucleotid is not the same as in the given sequence or can not be sequenced
                        if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("t" != choosenReadSeq[insertPlace]) and ("N" != choosenReadSeq[insertPlace]):

                            # inserting the mutations into the copied sequence
                            mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                            mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                            newSeq = mutableChoosenReadSeq.toseq()

                            newRecord = SeqRecord(newSeq, id=idName)

                            randomPostions.append(bins[bin]+insertPlace+1)
                            mutatedNucleotids.append(mutatedNucleotid)

                            inserted = 1

            # if there's no MDA amplification error, the choosen copy of the sequence gets saved into a new record
            else:
                newRecord = SeqRecord(choosenReadSeq, id =idName)

            # now we add the new record to the records vector
            records.append(newRecord)

        # after the amplification all the records get saved into one fasta file
        SeqIO.write(records, ("Data/MDASimulation/" + nametemplate), "fasta")

        # using the given read length, mean fragment size and standart deviation of fragment size art_illumina simulates
        # the illumina sequencing system  HiSeq 2500 with paired-end read simulation and user specified ART illumina sequence coverage
        outputfile = "results/ARTRecord_Tree" + str(snakemake.wildcards.treename) + "_Sample" + str(snakemake.wildcards.sample) + "_Bin" + str(bin) + "."
        shellCommand = "art_illumina -na -i " + ("Data/MDASimulation/" + nametemplate) + " -p -l " + str(readLength) + " -ss HS25 -f " + str(ARTcoverage) + " -m " + str(meanFragmentSize) + " -s " + str(standartDeviationOfFragmentSize) + " -o " + outputfile
        os.system(shellCommand)

        # make sure ART is finished before continuing
        time_to_wait = 1000
        time_counter = 0
        while not os.path.exists(outputfile + "1.fq"):
            time.sleep(1)
            time_counter += 1
            if time_counter > time_to_wait:break

        # the fasta file is only needed for the art_illumina simulation and is deleted to save space
        shellCommand = "rm Data/MDASimulation/" + nametemplate
        os.system(shellCommand)

        nameForFASTQFiles = "Chromosom " + str(chromosom) + " Allel " + str(snakemake.wildcards.sample)

        numberOfReads = 0
        for rec in SeqIO.parse(outputfile + "1.fq", "fastq"):
            numberOfReads += 1

        # choose with if a read gets lost during sequencing with some probability ARTDropReadProbability
        # 0 represents a read not being dropped and 1 a read being dropped
        droppedReads = random.choices(population = [0, 1], weights = [(1-ARTDropReadProbability), ARTDropReadProbability], k = numberOfReads)

        with open(outputfile + "1.fq") as in_handle:
            i = 0
            for rec in FastqGeneralIterator(in_handle):
                # if read i has a 0 in the droppedReads vector this means it hasn't been dropped and gets added to the final reads
                if (droppedReads[i] == 0):
                    ArtIlluminaRecords1.append(rec)
                i = i + 1

        # the art_illumina output file is not needed anymore and is deleted to save space
        shellCommand = "rm " + outputfile + "1.fq"
        os.system(shellCommand)

        with open(outputfile + "2.fq") as in_handle:
            i = 0
            for rec in FastqGeneralIterator(in_handle):
                # if read i has a 0 in the droppedReads vector this means it hasn't been dropped and gets added to the final reads
                if (droppedReads[i] == 0):
                    ArtIlluminaRecords2.append(rec)
                i = i + 1

        # the art_illumina output file is not needed anymore and is deleted to save space
        shellCommand = "rm " + outputfile + "2.fq"
        os.system(shellCommand)


# at the end all the choosen records from the resultReads lists are written into two resulting files
with open(snakemake.output[0], "w") as out_handle:
    for i in range(0, len(ArtIlluminaRecords1)):
        out_handle.write("@%s\n%s\n+\n%s\n" % (nameForFASTQFiles, ArtIlluminaRecords1[i][1], ArtIlluminaRecords1[i][2]))

with open(snakemake.output[1], "w") as out_handle:
    for i in range(0, len(ArtIlluminaRecords2)):
        out_handle.write("@%s\n%s\n+\n%s\n" % (nameForFASTQFiles, ArtIlluminaRecords2[i][1], ArtIlluminaRecords2[i][2]))

vcfReader_MetaData = vcf.Reader(filename= snakemake.params[13])
vcfWriter = vcf.Writer(open(vcfName, 'w'), vcfReader_MetaData)
chromosomId = int(next(iter(vcfReader_MetaData.contigs)))

for mutation in range(0, len(mutatedNucleotids)):

    position = int(randomPostions[mutation])
    insertMutation = mutatedNucleotids[mutation]
    record = vcf.model._Record(CHROM=chromosomId, POS=(position), ID='.',
                REF=vcf.model._Substitution(referenceGenome[position-1]),
                ALT=[vcf.model._Substitution(insertMutation)], QUAL='.', FILTER='PASS', INFO={},
                FORMAT=".", sample_indexes=[], samples=None)
    vcfWriter.write_record(record)
    vcfWriter.flush()

# close the writer
vcfWriter.close()
