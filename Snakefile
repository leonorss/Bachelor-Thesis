configfile: "config.json"
SAMPLES = range(0, (2*(config["generateBinaryTree"]["numberOfLeaves"])))

rule all:
    input:
        "results/SingleCellReadSimulator_{treename}.txt"

rule generateBinaryTree:
    output:
        "Data/Trees/binaryTree{treename}.txt"

    params: config["generateBinaryTree"]["numberOfLeaves"],
            config["generateBinaryTree"]["seed"]

    conda:
        "envs/environment1.yaml"

    shell:
        "python scripts/generateBinaryTree.py {params[0]} {output} {params[1]}"

rule plotBinaryTree:
    input:
        "Data/Trees/binaryTree{treename}.txt"

    output:
        "Data/Trees/Graph{treename}.pdf"

    params:
        config["generateBinaryTree"]["numberOfLeaves"]

    conda:
        "envs/environment2.yaml"

    script:
        "scripts/plotTree.py"

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
        "envs/environment3.yaml"

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
        "envs/environment4.yaml"

    script:
        "scripts/createMutationsFromTree.py"

rule sortMutations:
    input:
        "Data/Mutations/Mutations{treename}_{sample}.vcf"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf")

    conda:
        "envs/environment5.yaml"

    shell:
        "bcftools sort {input} > {output}"

rule zipMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz")

    conda:
        "envs/environment5.yaml"

    shell:
        "bgzip {input}"

rule indexMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz"

    output:
        temp("Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz.tbi")

    conda:
        "envs/environment5.yaml"

    shell:
        "tabix {input}"

rule insertMutations:
    input:
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz",
        "Data/Mutations/sortedMutations{treename}_{sample}.vcf.gz.tbi"

    output:
        temp("Data/Mutations/InsertedMutations{treename}_{sample}.fa")

    params:
        config["setReferenceGenomeAndMetaDataFrame"]["referenceGenome"]

    conda:
        "envs/environment5.yaml"

    shell:
        "cat {params} | bcftools consensus {input[0]} > {output}"

rule simulateReads:
    input:
        "Data/Mutations/InsertedMutations{treename}_{sample}.fa"

    output:
        "results/simulatedAmplificationAndSequencing_{treename}_Allel{sample}_1.fq",
        "results/simulatedAmplificationAndSequencing_{treename}_Allel{sample}_2.fq"

    params:
        config["simulateReads"]["seed"],
        config["simulateReads"]["meanBinSize"],
        config["simulateReads"]["standartDeviationOfBinSize"],
        config["simulateReads"]["meanMDACoverage"],
        config["simulateReads"]["standartDeviationOfMDACoverage"],
        config["simulateReads"]["readLength"],
        config["simulateReads"]["meanFragmentSize"],
        config["simulateReads"]["standartDeviationOfFragmentSize"],
        config["simulateReads"]["MDAamplificationErrorProbability"],
        config["setReferenceGenomeAndMetaDataFrame"]["metaDataFrame"],
        config["simulateReads"]["ARTcoverage"],
        config["simulateReads"]["ARTDropReadProbability"],
        config["simulateReads"]["zeroInflatationProbability"],
        config["setReferenceGenomeAndMetaDataFrame"]["metaDataFrame"],
        config["setReferenceGenomeAndMetaDataFrame"]["referenceGenome"]

    conda:
        "envs/environment6.yaml"

    script:
        "scripts/ARTSimulation.py"

rule SingleCellReadSimulator:
    input:
        expand("results/simulatedAmplificationAndSequencing_{{treename}}_Allel{sample}_1.fq", sample = SAMPLES),
        expand("results/simulatedAmplificationAndSequencing_{{treename}}_Allel{sample}_2.fq", sample = SAMPLES),
        "Data/Trees/Graph{treename}.pdf"

    output:
        "results/SingleCellReadSimulator_{treename}.txt"

    shell:
        "echo 'Simulation ran succesfully.' > {output}"

rule cleanAll:
    shell:
        "rm -r Data/Trees Data/Mutations Data/MDASimulation results"
