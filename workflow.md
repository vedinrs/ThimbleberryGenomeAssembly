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
cd /project/def-mtodesco/vschimma/packages
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
			-- scripts -> /project/def-mtodesco/vschimma/packages/juicer/CPU
			/fastq
                        	--hic_R1.fastq -> [path-to-fastq]/hic-reads-SRR30502875_1.fastq
                        	--hic_R2.fastq -> [path-to-fastq]/hic-reads-SRR30502875_2.fastq
			/references
 				-- thimbleberry.asm.hic.hap1.p_ctg.fa
                    		-- thimbleberry.asm.hic.hap1.p_ctg.fa.amb
                    		-- thimbleberry.asm.hic.hap1.p_ctg.fa.ann
                    		-- thimbleberry.asm.hic.hap1.p_ctg.fa.bwt
                    		-- thimbleberry.asm.hic.hap1.p_ctg.fa.pac
                    		-- thimbleberry.asm.hic.hap1.p_ctg.fa.sa
                	/restriction-sites -> ../restriction-sites/
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
export PATH=$PATH:/project/def-mtodesco/vschimma/packages/juicer/CPU

cd /project/def-mtodesco/vschimma/thimbleberry/juicer/

#run juicer
bash scripts/juicer.sh \
-g rp_hap1_ctg -s DpnII \
-S early \
-p restriction_sites/rp.chrom.sizes \
-y restriction_sites/rp_hap1_ctg_DpnII.txt \
-z references/thimbleberry.asm.hic.hap1.p_ctg.fa \
-t 32
-D /project/def-mtodesco/vschimma/packages/juicer
--assembly

# ---------------------------------------------------------------------
echo "Done Juicer Hi-C analysis.  Use yahs to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
```

Note: "-D /project/def-mtodesco/vschimma/packages/juicer" is incredibly important as you are directing juicer.sh to where all the other juicer commands are being stored.

## Prepare files for yahs

Use the following script to convert merge_dups.txt into a .bed file, allowing for easier usage of yahs.

```
#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --mem=12G
#SBATCH --account=def-mtodesco

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

cd ./juicer/aligned

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
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32
#SBATCH --account=def-mtodesco
#SBATCH --output=yahs.out
#SBATCH --error=yahs.err

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

export PATH=$PATH:/home/vschimma/packages/yahs/yahs

yahs -o yahs-outfiles/ -e GATC ./juicer/references/thimbleberry.asm.hic.hap1.p_ctg.fa ./juicer/aligned/merged_nodups_for_yahs.bed

# ---------------------------------------------------------------------
echo "Done yahs pipeline assembly. Created scaffolds from contigs and Hi-C heatmap."
# ---------------------------------------------------------------------

echo "Finished job at `date`"

```

### Create the hic and assembly files for Juicerbox Assembly Tools

In order to use the scaffold and hic map in Juicebox, prepare the appropriate files using this script:

```
#!/bin/bash

#SBATCH --time=3-0
#SBATCH --mem=256G
#SBATCH --account=def-mtodesco
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

module load StdEnv/2020 python/3.11.2 java/17.0.2 lastz/1.04.03

/home/vschimma/packages/yahs/juicer pre -a -o yahs-outfiles yahs.out.bin yahs.out_scaffolds_final.agp thimbleberry.asm.hic.hap1.p_ctg.fa >yahs-outfiles.log 2>&1

java -jar -Xmx240G /project/def-mtodesco/vschimma/thimbleberry/juicer/scripts/common/juicer_tools.1.9.9_jcuda.0.8.jar pre yahs-outfiles.txt yahs-outfiles.hic.part <(cat yahs-outfiles.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv yahs-outfiles.hic.part yahs-outfiles.hic)

echo "Finished job at `date`"
```

## Manually polish with Juicebox Assembly Tools

Create a directory called /juicebox-outfiles.



## Use BUSCO to assess genome assembly quality

Create a directory called /busco-outfiles
