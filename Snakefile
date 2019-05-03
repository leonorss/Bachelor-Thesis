configfile: "config.json"
SAMPLES = range(0, (2*(config["generateBinaryTree"]["numberOfLeaves"])))

rule all:
    input:
        "results/SingleCellReadSimulator_report{treename}.txt"

rule generateBinaryTree:
    output:
        "Data/Trees/binaryTree{treename}.txt"

    params: leaves = config["generateBinaryTree"]["numberOfLeaves"],
            seed = config["generateBinaryTree"]["seed"]

    shell:
        "python scripts/generateBinaryTree.py {params[0]} {output} {params[1]}"

rule createGermlineMutations:
    input:
        "Data/ReferenceData/ReferenceGenomeChromosome_22.fa",
        "Data/ReferenceData/MetaDataFrame.vcf"

    output:
        "Data/Mutations/GermlineMutations{treename}_Allel1.vcf",
        "Data/Mutations/GermlineMutations{treename}_Allel2.vcf"

    params:
        n_g_m = config["createGermlineMutations"]["numberOfGermlineMutations"],
        f_h_m = config["createGermlineMutations"]["fractionOfHeterozygousGermlineMutations"],
        seed = config["createGermlineMutations"]["seed"]

    script:
        "scripts/createGermlineMutations.py"

# rule insertGermlineMutations:
#     input:
#         genome = "Data/ReferenceData/ReferenceGenomeChromosome_18.fa",
#         mutationsAllel1 = "Data/Mutations/GermlineMutations{treename}Allel1.vcf",
#         mutationsAllel2 = "Data/Mutations/GermlineMutations{treename}Allel2.vcf"
#
#     output:
#         "Data/Mutations/InsertedGermlineMutations{treename}.fa"
#
#     shell:
#         "touch {output}"

#        "insertMutationsIntoGenome.sh -I {input.mutations} -g {input.genome} > {output}"

rule createMutations:
    input:
        "Data/Trees/binaryTree{treename}.txt",
        "Data/Mutations/GermlineMutations{treename}_Allel1.vcf",
        "Data/Mutations/GermlineMutations{treename}_Allel2.vcf",
        "Data/ReferenceData/ReferenceGenomeChromosome_22.fa",
        "Data/ReferenceData/MetaDataFrame.vcf"

    output:
        expand("Data/Mutations/Mutations{{treename}}_{sample}.vcf", sample = SAMPLES)

    params:
        config["createMutations"]["numberOfMutations"],
        config["generateBinaryTree"]["numberOfLeaves"],
        config["createMutations"]["seed"]

    script:
        "scripts/createMutationsFromTree.py"

rule sortMutations:
    input:
        "Data/Mutations/Mutations{treename}_{sample}.vcf"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf")

    shell:
        "bcftools sort {input} > {output}"

rule zipMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz")

    shell:
        "bgzip {input}"

rule indexMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz.tbi")

    shell:
        "tabix {input}"

rule insertMutations:
    input:
        "Data/ReferenceData/ReferenceGenomeChromosome_22.fa",
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz",
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz.tbi"

    output:
        "Data/Mutations/InsertedMutations{treename}_{sample}.fa"

    shell:
        "cat {input[0]} | bcftools consensus {input[1]} > {output}"

rule simulateReads:
    input:
        "Data/Mutations/InsertedMutations{treename}_{sample}.fa",
        "Data/ReferenceData/ReferenceGenomeChromosome_22.fa"

    output:
        "Data/MDASimulation/simulatedRead{treename}_{sample}txt"

    params:
        config["simulateReads"]["seed"],
        config["simulateReads"]["meanBinSize"],
        config["simulateReads"]["standardDeviationOfBinSize"],
        config["simulateReads"]["meanCoverage"],
        config["simulateReads"]["standardDeviationCoverage"],
        config["simulateReads"]["readLength"],
        config["simulateReads"]["meanFragmentSize"],
        config["simulateReads"]["standardDeviationOfFragmentSize"],
        config["simulateReads"]["MDAamplificationErrorProbability"]

    script:
        "scripts/splitGenomesIntoSegments.py"

rule convertFilesToSamtools:
    input:
        "Data/MDASimulation/simulatedRead{treename}_{sample}txt"

    output:
        "Data/MDASimulation/convertedToSamtools{treename}_{sample}.txt"

    shell:
        "touch {output}"



rule simulateMDAAmplificationErrors:
    input:
        "Data/MDASimulation/combinedReads{treename}_{sample}.txt"

    output:
        "results/amplificationErrors{treename}_{sample}.txt"

    shell:
        "touch {output}"

rule SingleCellReadSimulator:
    input:
        expand("results/amplificationErrors{{treename}}_{sample}.txt", sample = SAMPLES)

    output:
        "results/SingleCellReadSimulator_report{treename}.txt"

    shell:
        "touch {output}"
