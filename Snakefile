configfile: "config.json"
SAMPLES = range(0, (2*(config["generateBinaryTree"]["numberOfLeaves"])))

rule all:
    input:
        "results/SingleCellReadSimulator{treename}.txt"

rule generateBinaryTree:
    output:
        "Data/Trees/binaryTree{treename}.txt"

    params: config["generateBinaryTree"]["numberOfLeaves"],
            config["generateBinaryTree"]["seed"]

    conda:
        "envs/environment1.yaml"

    shell:
        "python scripts/generateBinaryTree.py {params[0]} {output} {params[1]}"

rule createGermlineMutations:
    output:
        "Data/Mutations/GermlineMutations{treename}_Allel1.vcf",
        "Data/Mutations/GermlineMutations{treename}_Allel2.vcf"

    params:
        config["createGermlineMutations"]["numberOfGermlineMutations"],
        config["createGermlineMutations"]["fractionOfHeterozygousGermlineMutations"],
        config["createGermlineMutations"]["seed"],
        config["setReferenceGenomeAndMetaDataFrame"]["referenceGenome"],
        config["setReferenceGenomeAndMetaDataFrame"]["metaDataFrame"]

    conda:
        "envs/environment2.yaml"

    script:
        "scripts/createGermlineMutations.py"

rule createMutations:
    input:
        "Data/Trees/binaryTree{treename}.txt",
        "Data/Mutations/GermlineMutations{treename}_Allel1.vcf",
        "Data/Mutations/GermlineMutations{treename}_Allel2.vcf"

    output:
        expand("Data/Mutations/Mutations{{treename}}_{sample}.vcf", sample = SAMPLES)

    params:
        config["createMutations"]["numberOfMutations"],
        config["generateBinaryTree"]["numberOfLeaves"],
        config["createMutations"]["seed"],
        config["setReferenceGenomeAndMetaDataFrame"]["referenceGenome"],
        config["setReferenceGenomeAndMetaDataFrame"]["metaDataFrame"]

    conda:
        "envs/environment3.yaml"

    script:
        "scripts/createMutationsFromTree.py"

rule sortMutations:
    input:
        "Data/Mutations/Mutations{treename}_{sample}.vcf"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf")

    conda:
        "envs/environment4.yaml"

    shell:
        "bcftools sort {input} > {output}"

rule zipMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz")

    conda:
        "envs/environment4.yaml"

    shell:
        "bgzip {input}"

rule indexMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz.tbi")

    conda:
        "envs/environment4.yaml"

    shell:
        "tabix {input}"

rule insertMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz",
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz.tbi"

    output:
        "Data/Mutations/InsertedMutations{treename}_{sample}.fa"

    params:
        config["setReferenceGenomeAndMetaDataFrame"]["referenceGenome"]

    conda:
        "envs/environment4.yaml"

    shell:
        "cat {params} | bcftools consensus {input[0]} > {output}"

rule simulateReads:
    input:
        "Data/Mutations/InsertedMutations{treename}_{sample}.fa"

    output:
        temp("results/amplificationErrors{treename}_{sample}_calculated.txt")

    params:
        config["simulateReads"]["seed"],
        config["simulateReads"]["meanBinSize"],
        config["simulateReads"]["varianceOfBinSize"],
        config["simulateReads"]["meanCoverage"],
        config["simulateReads"]["varianceOfCoverage"],
        config["simulateReads"]["readLength"],
        config["simulateReads"]["meanFragmentSize"],
        config["simulateReads"]["varianceOfFragmentSize"],
        config["simulateReads"]["MDAamplificationErrorProbability"]

    conda:
        "envs/environment5.yaml"

    script:
        "scripts/ARTSimulation.py"

rule SingleCellReadSimulator:
    input:
        expand("results/amplificationErrors{{treename}}_{sample}_calculated.txt", sample = SAMPLES)

    output:
        "results/SingleCellReadSimulator{treename}.txt"

    shell:
        "echo 'Simulation ran succesfully.' > {output}"
