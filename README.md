# *simMutAccumTE*

## Description
This program is a modified version of 'simpoolTE' from Adrion et al. (2019) and used by Variation in transposable element activity over spatial and temporal scales (Ho et al. 2020, In review) to estimate the false positive and false negative rates of detecting tranposable element mutations when using TEFLoN (Adrion et al. 2017).
SimMutAccumTE simulates a mutation accumulation (MA) experiment by inserting and/or deleting a number of transposable elements (TEs) into a diploid genome of a focal MA line while leaving the ancestral (ANC) genomes intact. After mutations are simulated in the given genome, pIRS (Hu et al. 2012) is used to generate paired-end reads for all lines.

## Requirements
- Python v3.7
- pIRS (https://github.com/galaxy001/pirs)
- seqtk (https://github.com/lh3/seqtk)

## Usage

```
usage: simMutAccumTE.py <required> [optional] 
  -wd WD          <Full path to working directory (include prefix for new folders)>
  -pirs PIRSPATH  <Full path to pIRS (including /pirs at the end)>
  -c REFCHROM     <Chromosome to simulate in fasta format (must contain only one sequence)>
  -b BED          <TE annotation bed file>
  -nmut NMUT      <Number of TE mutations to simulate in MA line>
  -tmut TYPEMUT   <Type of mutation to simulate [1, 2, 3, 4]>
  -nhet NHET      <Number of shared heterozygous TE sites to simulate>
  -ncl NCL        <Number of non-focal MA lines to simulate (clones of Ancestor)>
  -mnlen MNLEN    [Minimum length of TEs to insert and delete]
  -mxlen MXLEN    [Maximum length of TEs to insert and delete]
  -r RANDSEED     [Seed for random number generator]
  -snp SNPRATE    [Rate of heterozygosity for haploid genome]
  -x COV          [Coverage for pIRS to simulate]
  -rlen RLEN      [Read length for pIRS to simulate]
  -insz INSIZE    [Insert size for pIRS to simulate]

```

SimMutAccumTE inserts or delete TEs from the TE annotation file (<em>b</em>) onto the genome fasta file (<em>c</em>).

First, pIRS (Hu et al. 2012) is used to generate a diploid ancestor (ANC) genome by duplicated the given genome fasta file and adding SNPs at a specified rate (<em>snp</em>); the fasta files of each homolog of ANC is stored in pirsAnc/ as ANC.H1.snp.fa and ANC.H2.snp.fa with the list of simualted SNPS recored in ANC.H1.snp.lst and ANC.H2.snp.lst.
pIRS requires that fasta file does not contain N's. 

Subsequently <em>nmut</em> TE mutations of a specified type (<em>tmut</em>) are simulated. Type 1 are novel TE insertions (0->1 gain), simulated by insertion a heterozygous TE (i.e., on one homolog of the genome) into the MA line. Type 2 are TE deletion from an ancestrally homozygous TE site (2->1 loss), simulated by inserting a homozygous TE (i.e. insertion of both homologs) on ANC and a heterozygous TE on the MA line. Type 3 are TE insertion onto an ancestrally heterozygous TE site (1->2 gain), simualted by inserting a heterozygous TE on the ANC and then a homozygous TE on the MA line. Type 4 are TE deletions from an ancestrally heterozygous TE site (1->0 loss), simualted by inserting a heterozygous TE on the ANC only. In addition a <em>nhet</em> number of shared TE heterozygous sites can be added onto both ANC and MA lines. 

TEs to be inserted where chosen randomly from the annotation file with a minimum and maximum length given by <em>mnlen</em> and <em>mxlen</em>, respectively. Insertions sites are randomly chosen while ensuring that the inserted TE is separated from all other existing TEs by at least one read length (<em>rlen</em>). TE insertions sites are also flanked by a target site duplication with a mean length of 5 bp drawn from a Poisson distribution. Fasta files of the homologs containing the simulated mutations for ANC and MA are stored in simTE/ as ANC.H1.simTE.fa,  ANC.H2.simTE.fa,  MA.H1.simTE.fa and  MA.H2.simTE.fa

After all mutations have been simulated, pIRS (Hu et al. 2012) is used to generate paired-end reads with read length <em>rlen</em> and insert size <em>insz</em> at a diploid coverage of <em>x</em>. Reads for the ancestral and <em>ncl</em> non-focal descendent lines are simulated using ANC.H1.simTE.fa and ANC.H2.simTE.fa while those for the focal MA line are simulated from MA.H1.simTE.fa and MA.H2.simTE.fa. Reads for all lines are stored as fastq format in pirsReads/.

## References
- Adrion JR, Song MH, Schrider DR, Hahn MW, Schaack S. Genome-wide estimates of transposable element insertion and deletion rates in Drosophila melanogaster. Genome Biol Evol. 2017;9(5):1329-40.
- Adrion, JR, Begun, DJ, Hahn MW. Patterns of transposable element variation and clinality in Drosophila. Mol Biol. 2019;28:1523-36.
- Hu X, Yuan J, Shi Y, Lu J, Liu B, Li Z, Chen Y, Mu D, Zhang H, Yue Z, Bai F, Li H, Fan W. pIRS: Profile-based illumine pair-end reads simulator. Bioinformatics. 2012;28(11):1533-35.

## Citation
If you used simMutAccumTE in you work, please cite:
Ho EKH, Bellis ES, Calkins J, Adrion JR, Latta IV LC, Schaack S (2021) Engines of change: Transposable element mutation rates are high and variable within *Daphnia magna*. PLoS Genet 17(11): e1009827. [Link](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009827)
