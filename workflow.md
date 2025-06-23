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

Use FastQC to verify high quality reads. Create a new directory and cd into it to store FastQC outputs. Then, run FastQC on the newly downloaded fastq files in a job:

```
#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --account=
#SBATCH --out=fastqc-raw-reads.out

cd /project/def-mtodesco/vschimma/thimbleberry/

module load fastqc

fastqc -o ./fastqc-outfiles ./raw-data/*.fastq
```


## Draft genome assembly with hifiasm

Now, we will use hifiasm to change the HiFi reads into an assembly. 

This code is adapted from a [workflow](https://github.com/kaede0e/stinging_nettle_genome_assembly/blob/main/1_genome_assembly_with_PBHiFi_HiC/how_to_run_this_assembly_workflow.md) on genome assembly by Kaede H. 

If it is taking a long time to find a node, use 32 cpus instead of 48. Make sure to change it in both the SBATCH comment and in the hifiasm command.

Make a new directory called /hifiasm-outfiles as well.

```
#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=192000M
#SBATCH --account=def-mtodesco
#SBATCH --cpus-per-task=32
#SBATCH --output=hifiasm.out
#SBATCH --error=hifiasm.err

### Draft genome assembly with HiFiasm + Hi-C data ###

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

module load StdEnv/2020
module load hifiasm/0.19.5

#####################################
##### Variables / data ##############
#####################################

HiC_1=/project/def-mtodesco/vschimma/thimbleberry/raw-data/hic-reads-SRR30502875_1.fastq
HiC_2=/project/def-mtodesco/vschimma/thimbleberry/raw-data/hic-reads-SRR30502875_2.fastq
Hifi_fastq=/project/def-mtodesco/vschimma/thimbleberry/raw-data/raw-reads-SRR30533059.fastq

cd hifiasm-outfiles/

hifiasm -o thimbleberry.asm -t 32 --h1 $HiC_1 --h2 $HiC_2 $Hifi_fastq -s 0.32 --hom-cov 60

# ---------------------------------------------------------------------
echo "Done assembly with Hifiasm. Use 3D-DNA to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
```

The next few steps for the pipeline will be split up as both haplotypes are processed. Make a folder for each haplotype (hap1 and hap2 will be used in this pipeline). The process for haplotype 1 will be used in the scripts, but will be the same for haplotype 2.

## Prepare files for Juicer

### Use BWA to create references

Create a directory called /references in your haplotype folder. Convert the .gfa file of haploid 1 or 2 from the hifiasm output into a fasta file with awk (as per the [hifiasm FAQ](https://hifiasm.readthedocs.io/en/latest/faq.html)):

```
awk '/^S/{print ">"$2;print $3}' /hifiasm-outfiles/thimblerry.asm.hic.hap1.p_ctg.gfa > hap1/references/thimbleberry.asm.hic.hap1.p_ctg.fa
```

Next, cd into /references. Then, use BWA to index by running this job:

```
#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --account=def-mtodesco
#SBATCH --output=hap2/bwa-indexing.out

cd /project/def-mtodesco/vschimma/thimbleberry/hap2/references/

module load bwa

bwa index *.fa
```

### Git clone juicer

Cd into wherever you store packages, make a juicer directory, and git clone juicer:

```
cd /project/def-mtodesco/vschimma/packages
git clone https://github.com/aidenlab/juicer.git
```

### Generating site positions

Return back to the project directory. Create a new directory called /restriction_sites and cd into it. Then, copy generate_site_positions.py ([path-to-packages]/juicer/misc/generate_site_positions.py) into /restriction_sites. Modify it by adding a new entry in filenames:

```
  filenames = {
    'hg19': '/seq/references/Homo_sapiens_assembly19.fasta',
    'mm9' : '/seq/references/Mus_musculus_assembly9.fasta',
    'mm10': '/seq/references/Mus_musculus_assembly10.fasta',
    'hg18': '/seq/references/Homo_sapiens_assembly18.fasta',
    'rp_hap1_ctg': '../references/thimbleberry.asm.hic.hap1.p_ctg.fa', #MODIFY HERE!!! put your contig/assembly name and its path
  }
```

 Then, run the following job to generate the site positions:

```
#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --account=def-mtodesco
#SBATCH --output=hap2/generate-site-positions.out

cd /project/def-mtodesco/vschimma/thimbleberry/hap2/restriction_sites/

module load python

python generate_site_positions.py DpnII rp_hap2_ctg # use the name you made in generate_site_positions.py, not the file path
```

### Get the chromosome sizes

Find the sizes of all the chromosomes:

```
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --account=
#SBATCH --output=chromosome-sizes.out

cd /project/def-mtodesco/vschimma/thimbleberry/hap1/restriction_sites/

for i in $(ls *DpnII.txt); do
  name=$(echo $i | cut -d "." -f 1 )
  awk 'BEGIN{OFS="\t"}{print $1, $NF}'  $i > "$name"".chrom.sizes"
done
```

### Set up the fastq files

Make a directory called /juicer and inside make another directory called /fastq. After cd'ing into /fastq, make links to the hi-c reads in the /raw-reads directory using ln:
```
ln -s [path-to-read]/read_1.fastq read_R1.fastq
ln -s [path-to-read]/read_2.fastq read_R2.fastq
```

Make sure that the reads include _R1 and _R2 respectively, otherwise juicer.sh will not be able to interpret them correctly.

### Link to a juicer script

In the directory the job is being run in, make a link called scripts to the /CPU directory within the juicer package. 

```
ln -s /project/def-mtodesco/vschimma/packages/juicer/CPU
```

Next, cd into scripts/common/ and use set up a juicer_tools jar file:

```
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
```

Overall, the directory structure should look something like this:

```
[path-to-dir]/thimbleberry
		/juicer
			/fastq
                        	--hic_R1.fastq -> [path-to-fastq]/hic-reads-SRR30502875_1.fastq
                        	--hic_R2.fastq -> [path-to-fastq]/hic-reads-SRR30502875_2.fastq
			references -> ../references/
                	restriction-sites -> ../restriction-sites/
			scripts -> /project/def-mtodesco/vschimma/packages/juicer/CPU
```

## Ruin juicer

Use the following script to run juicer with many CPUs:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32
#SBATCH --account=def-mtodesco
#SBATCH --output=hap1/juicer.out
#SBATCH --error=hap1/juicer.err

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
export PATH=$PATH:/project/def-mtodesco/vschimma/packages/juicer/CPU

cd /project/def-mtodesco/vschimma/thimbleberry/hap1/juicer/

#run juicer
bash scripts/juicer.sh \
-g rp_hap1_ctg -s DpnII \
-S early \
-p restriction_sites/rp.chrom.sizes \
-y restriction_sites/rp_hap1_ctg_DpnII.txt \
-z references/thimbleberry.asm.hic.hap1.p_ctg.fa \
-t 32 \
-D $PWD \
--assembly

# ---------------------------------------------------------------------
echo "Done Juicer Hi-C analysis.  Use yahs to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
```

Note: "-D $PWD" is incredibly important as you are directing juicer.sh to where all the other juicer commands are being stored.

## Prepare files for yahs

Use the following script to convert merge_dups.txt into a .bed file, allowing for easier usage of yahs.

```
#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --mem=12G
#SBATCH --account=def-mtodesco
#SBATCH --output=hap1/make-bed.out

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

cd ./hap1/juicer/aligned

awk '
{
     if ($9 > 0 && $12 > 0) {
    # Calculate end position for read 1 using pos1 and cigar1
    cigar = $10
    pos_start1 = $3
    sum1 = 0
    while (match(cigar, /([0-9]+)([MDNX=])/, arr)) {
        sum1 += arr[1]
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    pos_end1 = pos_start1 + sum1  # End position for read 1

    # Calculate end position for read 2 using pos2 and cigar2
    cigar = $13
    pos_start2 = $7
    sum2 = 0
    while (match(cigar, /([0-9]+)([MDNX=])/, arr)) {
        sum2 += arr[1]
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    pos_end2 = pos_start2 + sum2  # End position for read 2

    # Print interleaved output for read 1 and read 2
    print $2, pos_start1, pos_end1, $15"/1", $9
    print $6, pos_start2, pos_end2, $16"/2", $12
}
}
' merged_nodups.txt > merged_nodups_for_yahs.bed

echo "Finished job at `date`"
```

## Scaffold assembly using yahs

Make a directory for the output from yahs called /yahs-outfiles. Then, run this job:

```
#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --account=def-mtodesco
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --output=hap1/yahs.out
#SBATCH --error=hap1/yahs.err


#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

# --- MODULES ---

module load StdEnv/2023 gcc/12.3 samtools/1.20

# ---------------

export PATH=$PATH:/home/vschimma/packages/yahs

CONTIG=.hap1/juicer/references/thimbleberry.asm.hic.hap1.p_ctg.fa

samtools faidx $CONTIG

yahs -o hap1/yahs-outfiles/yahs.out -e GATC $CONTIG ./hap1/juicer/aligned/merged_nodups_for_yahs.bed

# ---------------------------------------------------------------------
echo "Done yahs pipeline assembly. Created scaffolds from contigs and Hi-C heatmap."
# ---------------------------------------------------------------------

echo "Finished job at `date`"

```

### Create the hic and assembly files for Juicerbox Assembly Tools

Create a directory called juicebox-infiles. In order to use the scaffold and hic map in Juicebox, prepare the appropriate files using this script:

```
#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --account=def-mtodesco
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=hap1/hic-assembly.out

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

# --- MODULES ---

module load StdEnv/2020 python/3.11.2 java/17.0.2 lastz/1.04.03

# ---------------

cd hap1/juicebox-infiles/

/home/vschimma/packages/yahs/juicer pre -a -o JBAT ../yahs-outfiles/yahs.out.bin ../yahs-outfiles/_scaffolds_final.agp ../juicer/references/thimbleberry.asm.hic.hap1.p_ctg.fa.fai >JBAT.log 2>&1

(java -jar -Xmx240G /project/def-mtodesco/vschimma/thimbleberry/juicer/scripts/common/juicer_tools.1.9.9_jcuda.0.8.jar pre JBAT.txt JBAT.hic.part <(cat JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) \
&& (mv JBAT.hic.part JBAT.hic)

echo "Finished job at `date`"
```

## Manually polish with Juicebox Assembly Tools

Load the .hic file from juicebox-infiles/ and use the corresponding .assembly file in Juicebox Assembly Tools. Then, complete manual curation. The [Juicebox tutorial video](https://www.youtube.com/watch?v=Nj7RhQZHM18) may be helpful.

Once complete, create a directory called /juicebox-outfiles and save the result there.


## Use BUSCO to assess genome assembly quality

Create a directory called /busco-outfiles. Throughout the pipeline, it is good to use BUSCO to analyze the quality of the assembly. To do so, use the following script:

```

```

## Reorder genome based on raspberry

The next step is to align the chromosomes with the raspberry chromosomes to make annotation later on much easier. To start this process, we will complete some manual curation and use Mummer.

## Removing short sequences

At this point, there are some sequences not placed on any particular chromosome. It's important to remove these sequences to then align the rest of the data to the raspberry genome. To do this, use the following script:

```
#!/bin/bash
#SBATCH --account=def-mtodesco
#SBATCH --time=00:50:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --out=remove-short-seq.out
#SBATCH --err=remove-short-seq.err

# ---------------------------------------------------------------------
echo ""
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
echo ""
# ---------------------------------------------------------------------

# --- Modules ---

module load StdEnv/2020
module load seqkit/2.3.1
module load bioawk/1.0

# ---------------

cd ~/scratch/thimbleberry

# remove short scaffolds

seqkit seq -m 1000000 ./juicebox-outfiles/hap1-unaligned.fa > ./mummer-infiles/hap1.fasta
seqkit seq -m 1000000 ./juicebox-outfiles/hap2-unaligned.fa > ./mummer-infiles/hap2.fasta
seqkit seq -m 1000000 ./raspberry.fa > ./mummer-infiles/raspberry.fasta


# rename scaffolds

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  ./mummer-infiles/hap1.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">hap1-" ++i "-" length($seq) "\n" $seq}'  > ./mummer-infiles/hap1-names.fasta

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  ./mummer-infiles/hap2.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">hap2-" ++i "-" length($seq) "\n" $seq}'  > ./mummer-infiles/hap2-names.fasta

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  ./mummer-infiles/raspberry.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">raspberry-" ++i "-" length($seq) "\n" $seq}'  > ./mummer-infiles/raspberry-names.fasta

grep ">" ./mummer-infiles/*-names.fasta


# ---------------------------------------------------------------------
echo "Done prep for Mummer."
echo "Finished job at `date`"
# ---------------------------------------------------------------------
```

This script also renames the fasta entries for both haplotypes as well as for the raspberry sequences.

## Using Mummer

Use Mummer to produce an alignment plot using the following script:

```
#!/bin/bash
#SBATCH --account=def-mtodesco
#SBATCH --time=00:50:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --out=mummer.out
#SBATCH --err=mummer.err

# ---------------------------------------------------------------------
echo ""
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
echo ""
# ---------------------------------------------------------------------

# --- Modules ---

#module load StdEnv/2023
#module load python/3.10.2 minimap2/2.24 scipy-stack/2022a gcc/9.3.0 r-bundle-bioconductor/3.12 mummer/4.0.0beta2 seqtk
#module load seqkit

module load StdEnv/2023
module python minimap2 scipy-stack gcc r-bundle-bioconductor
module load seqtk
module load mummer


# ---------------

cd ~/scratch/thimbleberry

# run mummer

nucmer -t 40 -c 65 -l 20 -p ./mummer-outfiles/hap1 mummer-infiles/hap1.fasta  mummer-infiles/raspberry-names.fasta
nucmer -t 40 -c 65 -l 20 -p ./mummer-outfiles/hap2 mummer-infiles/hap2.fasta  mummer-infiles/raspberry-names.fasta

delta-filter -l 10000 -q -r ./mummer-outfiles/hap1.delta > ./mummer-outfiles/hap1-filter.delta
delta-filter -l 10000 -q -r ./mummer-outfiles/hap2.delta > ./mummer-outfiles/hap2-filter.delta

mummerplot ./mummer-outfiles/hap1-filter.delta -R ./mummer-infiles/hap1.fasta -Q ./mummer-infiles/raspberry-names.fasta --png -p ./mummer-outfiles/mum-plot-hap1
mummerplot ./mummer-outfiles/hap2-filter.delta -R ./mummer-infiles/hap2.fasta -Q ./mummer-infiles/raspberry-names.fasta --png -p ./mummer-outfiles/mum-plot-hap2

# ---------------------------------------------------------------------
echo "Done Mummer."
echo "Finished job at `date`"
# ---------------------------------------------------------------------
```

Transfer this file to your local machine to view it. It may look something like this:

![mummer-plot](https://github.com/vedinrs/ThimbleberryGenomeAssembly/blob/main/mum-plot-hap1.png "Mummer Alignment Plot Example")

## Reordering chromosomes manually

The following commands were used to manually reorder the chromosomes within the fasta file.

```
seqkit split -i -O reordered/ ./mummer-infiles/hap1-names.fasta

seqkit seq -t DNA --reverse hap1-names.part_hap1-3-34435378.fasta > hap1-3-reverse.fasta
```

## Mask duplicated sequences using Red

Using Red, mask duplicated sequences using the folliwing script:

```

```

## Genome annotation with EDTA

To use EDTA, it is easiest to set up a Conda environment.
