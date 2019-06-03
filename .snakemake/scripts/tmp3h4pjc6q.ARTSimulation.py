
######## Snakemake header ########
import sys; sys.path.extend(["/usr/local/lib/python3.7/site-packages", "/Users/leonorschubert/Desktop/Aramis/ETH/Bachelor Thesis/Code/scripts"]); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X(\x00\x00\x00Data/Mutations/InsertedMutations10_13.faq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0cX$\x00\x00\x00results/amplificationErrors10_13.txtq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12(M\xdb\x1dK\xc8G?\xe6ffffffK\xfaG?\xe6ffffffKdK\xfaG?\xe6ffffffG>z\xd7\xf2\x9a\xbc\xafHe}q\x13h\x08}q\x14sbX\t\x00\x00\x00wildcardsq\x15csnakemake.io\nWildcards\nq\x16)\x81q\x17(X\x02\x00\x00\x0010q\x18X\x02\x00\x00\x0013q\x19e}q\x1a(h\x08}q\x1b(X\x08\x00\x00\x00treenameq\x1cK\x00N\x86q\x1dX\x06\x00\x00\x00sampleq\x1eK\x01N\x86q\x1fuX\x08\x00\x00\x00treenameq h\x18X\x06\x00\x00\x00sampleq!h\x19ubX\x07\x00\x00\x00threadsq"K\x01X\t\x00\x00\x00resourcesq#csnakemake.io\nResources\nq$)\x81q%(K\x01K\x01e}q&(h\x08}q\'(X\x06\x00\x00\x00_coresq(K\x00N\x86q)X\x06\x00\x00\x00_nodesq*K\x01N\x86q+uh(K\x01h*K\x01ubX\x03\x00\x00\x00logq,csnakemake.io\nLog\nq-)\x81q.}q/h\x08}q0sbX\x06\x00\x00\x00configq1}q2(X\x12\x00\x00\x00generateBinaryTreeq3}q4(X\x0e\x00\x00\x00numberOfLeavesq5K\nX\x04\x00\x00\x00seedq6MY\x01uX\x17\x00\x00\x00createGermlineMutationsq7}q8(X\x19\x00\x00\x00numberOfGermlineMutationsq9K\x14X\'\x00\x00\x00fractionOfHeterozygousGermlineMutationsq:G?\xd0\x00\x00\x00\x00\x00\x00h6K8uX\x0f\x00\x00\x00createMutationsq;}q<(X\x11\x00\x00\x00numberOfMutationsq=K\x01h6M\xea\x01uX\r\x00\x00\x00simulateReadsq>}q?(X\x0b\x00\x00\x00meanBinSizeq@K\xc8X\x1a\x00\x00\x00standardDeviationOfBinSizeqAG?\xe6ffffffh6M\xdb\x1dX\x0c\x00\x00\x00meanCoverageqBK\xfaX\x19\x00\x00\x00standardDeviationCoverageqCG?\xe6ffffffX\n\x00\x00\x00readLengthqDKdX\x10\x00\x00\x00meanFragmentSizeqEK\xfaX\x1f\x00\x00\x00standardDeviationOfFragmentSizeqFG?\xe6ffffffX \x00\x00\x00MDAamplificationErrorProbabilityqGG>z\xd7\xf2\x9a\xbc\xafHuuX\x04\x00\x00\x00ruleqHX\r\x00\x00\x00simulateReadsqIX\x0f\x00\x00\x00bench_iterationqJNX\t\x00\x00\x00scriptdirqKXE\x00\x00\x00/Users/leonorschubert/Desktop/Aramis/ETH/Bachelor Thesis/Code/scriptsqLub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/leonorschubert/Desktop/Aramis/ETH/Bachelor Thesis/Code/scripts/ARTSimulation.py';
######## Original script #########
import numpy as np
import random
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord

# read in seed from config file for reproducibilaty
configSeed = snakemake.params[0] + int(snakemake.wildcards.sample)
random.seed(configSeed)

# read in mean and variance of segment length from config file
meanBinSizeLength = snakemake.params[1]
standardDeviationOfBinSize = snakemake.params[2]

# read in genome sequence
genome = SeqIO.read(snakemake.input[0], "fasta")
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
while (position < (genomeLength - 2*meanBinSizeLength)):
    length = np.random.normal(meanBinSizeLength, standardDeviationOfBinSize)
    length = int(np.round(length))

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

# now we create the according fasta files
for bin in range(0, (len(bins)-1)):
    seqRecord = SeqRecord((genome[(bins[bin]):(bins[bin+1])].seq), id=(str(bin) + "." + str(0)))
    nametemplate = "Bin_Tree" + snakemake.wildcards.treename + "_Sample" + snakemake.wildcards.sample + "_Bin" + str(bin) + ".fa"
##    name = nametemplate + "0.fa"
##    SeqIO.write(seq_Record, name, "fasta")

    records = [seqRecord]
    readCoverage = np.random.normal(meanCoverage, standardDeviationCoverage)
    readCoverage = int(np.round(readCoverage))

    binLength = bins[bin+1] - bins[bin]

    for amplification in range(1, readCoverage):
        choosenReadNumber = random.randrange(0, (amplification))
        amplificationError = random.choices(population = [0, 1], weights = [(1-MDAamplificationErrorProbability), MDAamplificationErrorProbability], k = 1)

        choosenReadSeq = (records[choosenReadNumber]).seq
        idName = str(bin) + "." + str(amplification)

        if (amplificationError == 1):
            inserted = 0

            # checking that we didn't create a mutation that din't actually change the reference genome
            while(inserted == 0):
                # choosing random Nucleotid and random place to insert
                newMutation = random.randrange(1, 5)
                insertPlace = random.randrange(0, binLength)

                # case nucleotid="A"
                if newMutation==1:
                    mutatedNucleotid = "A"

                    # making sure the nucleotid is not the same as in the reference genome
                    if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("a" != choosenReadSeq[insertPlace]):

                    # inserting the mutations into the new file
                        mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                        mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                        newSeq = mutableChoosenReadSeq.toseq()

                        newRecord = SeqRecord(newSeq, id=idName)

                        inserted = 1




                # case nucleotid="C"
                elif newMutation==2:
                    mutatedNucleotid = "C"

                    # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
                    if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("c" != choosenReadSeq[insertPlace]):

                    # inserting the mutations into the new file
                        mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                        mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                        newSeq = mutableChoosenReadSeq.toseq()

                        newRecord = SeqRecord(newSeq, id=idName)

                        inserted = 1



                # case nucleotid="G"
                elif newMutation==3:
                    mutatedNucleotid = "G"

                    # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
                    if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("g" != choosenReadSeq[insertPlace]):

                    # inserting the mutations into the new file
                        mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                        mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                        newSeq = mutableChoosenReadSeq.toseq()

                        newRecord = SeqRecord(newSeq, id=idName)

                        inserted = 1


                # case nucleotid="T"
                else:
                    mutatedNucleotid = "T"

                    # making sure the nucleotid is not the same as in the reference genome or that we mutated the same nucleotide twize
                    if (mutatedNucleotid != choosenReadSeq[insertPlace]) and ("t" != choosenReadSeq[insertPlace]):

                    # inserting the mutations into the new file
                        mutableChoosenReadSeq = (choosenReadSeq).tomutable()
                        mutableChoosenReadSeq[insertPlace] = mutatedNucleotid
                        newSeq = mutableChoosenReadSeq.toseq()

                        newRecord = SeqRecord(newSeq, id=idName)

                        inserted = 1

        else:
            newRecord = SeqRecord(choosenReadSeq, id =idName)

        records.append(newRecord)

    SeqIO.write(records, nametemplate, "fasta")

    outputfile = "ARTRecord_Tree" + str(snakemake.wildcards.treename) + "_Sample" + str(snakemake.wildcards.sample) + "_Bin" + str(bin)
    shellCommand = "art_illumina -na -i " + nametemplate + " -p -l " + str(readLength) + " -ss HS25 -f 3 -m " + str(meanFragmentSize) + " -s " + str(standardDeviationOfFragmentSize) + " -o " + outputfile
    os.system(shellCommand)

shellCommand = "touch " + snakemake.output[0]
os.system(shellCommand)
