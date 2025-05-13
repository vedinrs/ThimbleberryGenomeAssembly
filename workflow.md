## Download from NCBI

To download the thimbleberry genome from NCBI, go to the SRAs from the primary haplotype ([PacBio Reads](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR30533059&display=download)) 
([Hi-C Reads](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR30502875&display=download)), and use the SRA toolkit. 

The SRR codes for this project are SRR30533059 and SRR30502875 for PacBio and Hi-C reads. To download the .sra files, use

```
prefetch SRRXXXXXXXX
```
and to convert them to fastq use
```
fasterq-dump --split-reads SRRXXXXXXXX.sra
```
while in the same directory as the files you just downloaded.


## Quality assurance

Use FastQC to verify high quality reads. Create a new directory and cd into it to store FastQC outputs. Then, run FastQC on the newly downloaded fastq files:

```
module load fastqc

fastqc [path-to-dir]/raw-data/*.fastq
```


## Draft genome assembly with hifiasm

Now, we will use hifiasm to change the HiFi reads into an assembly. 

This code is adapted from a [workflow](https://github.com/kaede0e/stinging_nettle_genome_assembly/blob/main/1_genome_assembly_with_PBHiFi_HiC/how_to_run_this_assembly_workflow.md) on genome assembly by Kaede H. 

```

```
