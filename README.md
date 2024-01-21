# EV - RNA Pipeline
We present a novel RNA sequencing (whole transcriptomics) analysis pipeline that quantifies gene/transcript counts from raw RNA sequencing reads based on alignment to the target organism's genome. This pipeline leverages on Linux, python and R packages to identify transcript packaging. Our pipeline integrates packages from the aforementioned tools for quality control, data pre-processing, differential gene analysis, and provides a comprehensive framework for RNA packaging studies. We demonstrate this pipeline’s performance and applicability to datasets from different organisms and biofluids.

## Nextflow
Nextflow is a workflow management system that allows for writing scalable and reproducible pipelines for automating analysis. It simplifies the execution of several complex tasks across computing platforms and environments. Singularity is a container manager that enables the running of applications in isolated and portable environments. It helps avoid conflicts of dependencies and ensures reproducibility results. Anaconda is a software distribution that provides a large repertoire of tools and libraries for data analysis, including sequencing-specific packages. It also allows the creation and management of virtual environments for different projects. Together, these software are deployed to provide a flexible, robust and efficient framework for RNA sequencing analysis.

# Installation
## Anaconda
We suggest using anaconda, which is a distribution of packages built for data science. Anaconda comes with conda, a package, and an environment manager that can be used to create environments for isolating projects that use different versions of Python and it’s packages.

1. Download the Anaconda installer for Linux from the [official website](https://docs.anaconda.com/anaconda/install/index.html).

2. Open a terminal window and navigate to the directory where you downloaded the installer.

3. Run the following command to add executable permission to the installer:
    ```bash
    chmod +x Anaconda3-YOUR_VERSION-Linux-x86_64.sh
    ```

4. Run the installer by running the following command:
    ```bash
    ./Anaconda3-YOUR_VERSION-Linux-x86_64.sh
    ```

5. Follow the prompts on the installer screens.

6. Once installation is complete, you can start using Anaconda by opening a new terminal window.

## Nextflow
Before installing Nextflow, you will need to make sure that **Java version 11** or greater is installed on your machine. You can check your Java version by running the following command in your terminal window:
```bash
java -version
```

1. Download the Nextflow executable by copying and pasting one of the following commands in your terminal window:
    ```bash
    wget -qO- https://get.nextflow.io | bash
    ```
    or
    ```bash
    curl get.nextflow.io | bash
    ```

2. Make the binary executable on your system by running:
    ```bash
    chmod +x nextflow
    ```

3. Optionally, move the nextflow file to a directory accessible by your $PATH:
    ```bash
    sudo mv nextflow /usr/local/bin/
    ```

### Nextflow with Conda

4. If you want to install Nextflow through Anaconda, you can do so by running:
    ```bash
    conda create -n nextflow -c bioconda nextflow
    ```

## Singularity
Singularity is a free and open-source computer program that performs operating-system-level virtualisation also known as containerization. One of the main uses of Singularity is to bring containers and reproducibility to scientific computing and the high-performance computing (HPC) world. Singularity is a container framework designed to run scientific applications on HPC- based resources.

Singularity is chosen by most HPCs as their primary container software because it allows users to pack an application/workflow/pipeline and all of its dependencies into a single image (file). This file can be easily transported between different HPC systems. Furthermore, Singularity assumes that the user does not have root privileges on the host operating system (OS), making it more secure than other containerization technologies.

1. Install dependencies:
    ```bash
    sudo apt-get update
    sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev wget pkg-config git cryptsetup-bin
    ```

2. Install Go:
    - Download the Go binary for Linux by going to the site [go.dev](https://go.dev/) and then clicking on Download Go.

    - Extract the Golang binaries tarball using the tar command to a directory of your choice:
        ```bash
        tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz
        ```
    - Extract the Golang binaries tarball using the tar command to a directory of your choice:
        ```bash
        export PATH=$PATH:/usr/local/go/bin
        ```
    - Extract the Golang binaries tarball using the tar command to a directory of your choice:
        ```bash
        go version
        ```

3. Download Singularity from a release:
    ```bash
    VERSION=3.8.0
    wget https://github.com/hpcng/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
    tar -xzf singularity-${VERSION}.tar.gz
    cd singularity-${VERSION}
    ```

4. Compile Singularity:
    ```bash
    ./configure --prefix=/usr/local
    make
    sudo make install
    ```

5. Verify that Singularity is installed correctly:
    ```bash
    singularity --version
    ```

### Singularity on the CAIR HPC

Alternatively if you are using CAIR you can enter the following commands to activate singularity:

1. load the GO module:
    ```bash
    module load application/go/1.14.2
    ```

2. load the Singularity module:
    ```bash
    module load application/singularity/3.5.3 
    ```
    *NOTE: The only version of singularity available on CAIR is 3.5.3, this should not cause issues with the execution of the pipeline.

# Running the Pipeline
## Minimal

```bash
conda activate nextflow
nextflow run /research/project/shared/benoukraf_lab/nextflow/microbiome_pipeline/main.nf -profile chia,ont,conda --input fastq/00811_8.fastq --target_index "/research/project/shared/benoukraf_lab/pathoscope/metascope/index/refseq_all/ont/*-ont.mmi" --filter_index /research/project/shared/benoukraf_lab/pathoscope/metascope/index/human/ont/hg38-ont.mmi
```

## Complex

```bash
conda activate nextflow
nextflow run Matthew-Dyer792/microbiome_pipeline -profile chia,ont,conda --conda_cacheDir /research/project/shared/benoukraf_lab/matthew/.conda_cacheDir --input fastq/00811_8.fastq --sequence_summary "summary/*.txt" --target_index "/research/project/shared/benoukraf_lab/pathoscope/metascope/index/refseq_all/ont/*-ont.mmi" --filter_index /research/project/shared/benoukraf_lab/pathoscope/metascope/index/human/ont/hg38-ont.mmi --ncbi_key t28545asd34jjhgjg2342
```

# Pipeline Options

### Pipeline location:
This is a path that leads to the folder containing the Microbiome Pipeline. It can be either a local path i.e. `"/research/project/shared/benoukraf_lab/nextflow/microbiome_pipeline"` or a github repository `"Matthew-Dyer792/microbiome_pipeline"`.
```bash
run Matthew-Dyer792/microbiome_pipeline
```

### Pipeline profile(s):
These are preset configurations that allow for standardize execution of the pipeline. For example there is a profile for conda and singularity, chia vs. local execution, and for oxford nanopore vs. illumina reads.
```bash
-profile illumina
```
or
```bash
-profile ont,conda,chia
```

### Pipeline workflow:
There are two workflow options with this pipeline, the primary is to run the metascope package on illumina reads. The secondary workflow recapitulates the metascope package but for oxford nanopore long reads. They can be selected with either illumina or ont. Note: Using the ont profile will set the workflow value to ont by default.
```bash
--workflow illumina
```

### Pipeline conda cache directory:
By default nextflow will store the conda files in the "work" directory. This file directory is usually discarded upon successful completion of the pipeline. Therefore if you want the conda directories to persist a external directory must be supplied. Note: If you intend to run this pipeline multiple times this will save the time and hard drive space of duplicate installations.
```bash
--conda_cacheDir /research/project/conda
```

### Input:
The primary input of the pipeline is fastq files. They can be provided as a regex "/data/big/\*.fq" or for paired end reads "data/big/file_{1,2}.fq". Note: the path must be encapsulated by "" in order to be properly interpreted.
```bash
--input "/research/project/fastq/*_{1,2}.fq"
```

### Single end
By default the pipeline assumes that paired end reads are provided. If single end reads are used set the flag to "true". Note: Using the ont profile will set single end to true by default.
```bash
--single_end true
```

### Sequence summary
The sequencing_summary.txt files created along with the oxford nanopore sequencing files. These are required for pycoqc. They can be provided as a regex "/data/big/\*.txt". Note: The path must be encapsulated by "" in order to be properly interpreted.
```bash
--sequence_summary "/research/project/summary/*.txt"
```

### Target index:
The genomes of microbes you aim to discover in your sample indexed by your aligner of choice (bowtie2, minimap2, bwa...). They can be provided as a regex "/data/big/\*.fq". Note: the path must be encapsulated by "" in order to be properly interpreted.
```bash
--target_index "/research/project/target_genomes/*-ont.mmi"
```

### Filter index:
The genomes of the host organism (or confounding organisms) you aim to remove from your sample indexed by your aligner of choice (bowtie2, minimap2, bwa...). They can be provided as a regex "/data/big/\*.mmi". Note: the path must be encapsulated by "" in order to be properly interpreted.
```bash
--filter_index "/research/project/filter_genomes/*-ont.mmi"
```

### NCBI key:
NCBI Entrez API key. optional. Due to the high number of requests made to NCBI, the ID function will be less prone to errors if you obtain an NCBI key. You may enter the string from your account.
```bash
--ncbi_key t28545asd34jjhgjg2342
```

### Pipeline resume:
A built in nextflow feature to allow resumption of the pipeline using a cached version of the previously completed steps. Note: "work" directory created by the pipeline must be present for this to function.
```bash
-resume
```