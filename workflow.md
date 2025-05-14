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

```
#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --mem=192000M
#SBATCH --account=
#SBATCH --cpus-per-task=48
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

HiC_1=/project/def-mtodesco/vschimma/thimbleberry/raw-data/
HiC_2=/project/def-mtodesco/vschimma/thimbleberry/raw-data/
Hifi_fastq=/project/def-mtodesco/vschimma/thimbleberry/raw-data/

hifiasm -o thimbleberry.asm -t 48 --h1 $HiC_1 --h2 $HiC_2 $Hifi_fastq -s 0.4 --hom-cov 128

# ---------------------------------------------------------------------
echo "Done assembly with Hifiasm. Use 3D-DNA to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
```

Make a new directory called /hifiasm-outfiles and place all the output files from hifiasm into it.

## Prepare files for Juicer

### Use BWA to create references

Create a directory called /references. Convert the .gfa file of haploid 1 from the hifiasm output into a fasta file with awk (as per the [hifiasm FAQ](https://hifiasm.readthedocs.io/en/latest/faq.html)):

```
awk '/^S/{print ">"$2;print $3}' /hifiasm-outfiles/thimblerry.asm.hic.hap1.p_ctg.gfa > references/thimbleberry.asm.hic.hap1.p_ctg.fa
```

Next, cd into /references. Then, use BWA to index by running this job:

```
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --account=

cd /project/def-mtodesco/vschimma/thimbleberry/references/

module load bwa

bwa index *.fa
```

### Git clone juicer

Cd into wherever you store packages, make a juicer directory, and git clone juicer:

```
cd /home/vschimma/packages/
mkdir juicer/
cd juicer/
git clone https://github.com/aidenlab/juicer.git
```

### Generating site positions

Return back to the project directory. Create a new directory called /restriction-sites and cd into it. Then, copy generate_site_positions.py ([path-to-packages]/juicer/misc/generate_site_positions.py) into /restriction-sites. Modify it by adding a new entry in filenames:

```
  filenames = {
    'hg19': '/seq/references/Homo_sapiens_assembly19.fasta',
    'mm9' : '/seq/references/Mus_musculus_assembly9.fasta',
    'mm10': '/seq/references/Mus_musculus_assembly10.fasta',
    'hg18': '/seq/references/Homo_sapiens_assembly18.fasta',
    'rp_hap1_ctg': '../references/thimbleberry.asm.hic.hap1.p_ctg.fa', #MODIFY HERE: put your contig/assembly and its path
  }
```

 Then, run the following job to generate the site positions:

```
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --account=

module load python

python generate_site_positions.py DpnII rp_hap1_ctg # use the name you made in generate_site_positions.py, not the file path
```

### Get the chromosome sizes

Find the sizes of all the chromosomes:

```
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --account=

for i in $(ls *DpnII.txt); do
  name=$(echo $i | cut -d "." -f 1 )
  awk 'BEGIN{OFS="\t"}{print $1, $NF}'  $i > "$name"".chrom.sizes"
done
```

## Ruin juicer

Use the following script to run juicer with many CPUs:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=510G
#SBATCH --cpus-per-task=32
#SBATCH --account=def-mtodesco
#SBATCH --output=juicer.out
#SBATCH --error=juicer.err

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
export PATH=$PATH:/home/~bin/juicer/CPU

PROJECT=/project/def-mtodesco/vschimma/thimbleberry/
cd $PROJECT/juicer-outfiles/

#run juicer
bash /home/vschimma/packages/juicer/CPU/juicer.sh \
-g rp_hap1_ctg -s DpnII \ #genome ID defined by -z command later, site
-S early \ # stage
-p $PROJECT/restriction_sites/thimbleberry.chrom.sizes \ # chromosome sizes
-y $PROJECT/restriction_sites/thimbleberry.asm.hic.hap1.p_ctg_DpnII.txt -z $PROJECT/references/thimbleberry.asm.hic.hap1.p_ctg.fa \ # restriction sites, reference genome
-t 28 # threads

# ---------------------------------------------------------------------
echo "Done Juicer Hi-C analysis.  Use 3D-DNA to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
```
