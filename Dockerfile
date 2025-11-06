# Use an existing image as the base
FROM ubuntu:latest

# Install necessary dependencies
RUN apt-get update && apt-get install -y \
    trimmomatic \
    fastqc \
    STAR \
    subread \
    wget \
    unzip \
    openjdk-8-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Download and install Nextflow
RUN wget -qO- https://get.nextflow.io | bash

# # Copy the Nextflow script into the container
# COPY rnaseq_script.nf /root/pipeline.nf

# # Copy the necessary reference files into the container
# COPY genomeIndex/ /root/genomeIndex/
# COPY annotation.gtf /root/annotation.gtf
# COPY adapters.fa /root/adapters.fa

# Set the working directory
WORKDIR /root

# Run the Nextflow script
# CMD ["nextflow", "run", "rna_script.nf"]