# nextPARS, a novel Illumina-based implementation of in-vitro parallel probing of RNA structures.

Here you will find the scripts necessary to produce the scores described in our paper from fastq files obtained during the experiment.

### Install Prerequisites
First install git:
```bash
sudo apt-get update
sudo apt-get install git-all
```

Then clone this repository

```bash
git clone https://github.com/Gabaldonlab/nextPARS.git
```

Now, ensure the necessary python packages are installed, and can be found in the `$PYTHONPATH` environment variable by running the script packages_for_nextPARS.sh in the nextPARS directory.

```bash
cd nextPARS/conf
chmod 775 packages_for_nextPARS.sh
./packages_for_nextPARS.sh
```


### Convert fastq to tab
In order to go from the fastq outputs of the nextPARS experiments to a format that allows us to calculate scores, first map the reads in the fastq files to a reference using the program of your choice.
Once you have obtained a bam file, use [PARSParser_0.67.b.jar](https://github.com/Gabaldonlab/nextPARS/tree/master/bin/PARSParser_0.67.b.jar).
This program counts the number of reads beginning at each position (which indicates a cut site for the enzyme in the file name) and outputs it in .tab format (count values for each position are separated by semi-colons).

Example usage:
```bash
java -jar PARSParser_0.67.b.jar -a bamFile -b bedFile -out outFile -q 20 -m 5
```

where the required arguments are:
  * -a gives the bam file of interest
  * -b is the bed file for the reference
  * -out is the name given to the output file in .tab format

Also accepts arguments: 
  * -q for minimum mapping quality for reads to be included [default = 0]
  * -m for minimum average counts per position for a given transcript [default = 5.0]



### Sample Data
There are sample data files found in the folder [nextPARS/data](https://github.com/Gabaldonlab/nextPARS/tree/master/data), as well as the necessary fasta files in [nextPARS/data/SEQS/PROBES](https://github.com/Gabaldonlab/nextPARS/tree/master/data/SEQS/PROBES), and the reference structures obtained from PDB in [nextPARS/data/STRUCTURES/REFERENCE_STRUCTURES](https://github.com/Gabaldonlab/nextPARS/tree/master/data/STRUCTURES/REFERENCE_STRUCTURES)
There are also 2 folders of sample output files from the PARSParser_0.67.b.jar program that can be used as further examples of the nextPARS score calculations described below. These folders are found in [nextPARS/data/PARSParser_outputs](https://github.com/Gabaldonlab/nextPARS/tree/master/data/PARSparser_outputs).
NOTE: these are randomly generated sequences with random enzyme values, so they are just to be used as examples for the usage of the scripts, good results should not be expected with these.


### nextPARS Scores
To obtain the scores from nextPARS experiments, use the script [get_combined_score.py](https://github.com/Gabaldonlab/nextPARS/tree/master/bin/get_combined_score.py). Sample data for the 5 PDB control structures can be found in the folder nextPARS/data/

There are a number of different command line options in the script, many of which were experimental or exploratory and are not relevant here. The useful ones in this context are the following:
  * Use the -i option [REQUIRED] to indicate the molecule for which you want scores (all available data files will be included in the calculations -- molecule name must match that in the data file names)
  * Use the -inDir option to indicate the directory containing the .tab files with read counts for each V1 and S1 enzyme cuts
  * Use the -f option to indicate the path to the fasta file for the input molecule

  * Use the -s option to produce an output Structure Preference Profile (SPP) file. Values for each position are separated by semi-colons. Here 0 = paired position, 1 = unpaired position, and NA = position with a score too low to determine its configuration.
  * Use the -o option to output the calculated scores, again with values for each position separated by semi-colons.
  * Use the --nP_only option to output the calculated nextPARS scores before incorporating the RNN classifier, again with values for each position separated by semi-colons.
  * Use the option {-V nextPARS} to produce an output with the scores that is compatible with the structure visualization program [VARNA](http://varna.lri.fr/)<sup>1</sup>
  * Use the option {-V spp} to produce an output with the SPP values that is compatible with VARNA.
  * Use the -t option to change the threshold value for scores when determining SPP values [default = 0.8, or -0.8 for negative scores]
  * Use the -c option to change the percentile cap for raw values at the beginning of calculations [default = 95]
  * Use the -v option to print some statistics in the case that there is a reference CT file available ( as with the example molecules, found in [nextPARS/data/STRUCTURES/REFERENCE_STRUCTURES](https://github.com/Gabaldonlab/nextPARS/tree/master/data/STRUCTURES/REFERENCE_STRUCTURES) ). If not, will still print nextPARS scores and info about the enzyme .tab files included in the calculations.

Example usage:
```bash
# to produce an SPP file for the molecule TETp4p6
python get_combined_score.py -i TETp4p6 -s
# to produce a Varna-compatible output with the nextPARS scores for one of the 
# randomly generated example molecules
python get_combined_score.py -i test_37 -inDir nextPARS/data/PARSParser_outputs/test1 \
  -f nextPARS/data/PARSParser_outputs/test1/test1.fasta -V nextPARS
```



### RNN classifier (already incorporated into the nextPARS scores above)
To run the RNN classifier separately, using a different experimental score input (in .tab format), it can be run like so:

```bash
python predict2.py -f molecule.fasta -p scoreFile.tab -o output.tab
```

Where the command line options are as follows:
  * the -f option [REQUIRED] is the input fasta file
  * the -p option [REQUIRED] is the input Score tab file
  * the -o option [REQUIRED] is the final Score tab output file.
  * the -w1 option is the weight for the RNN score. [default = 0.5]
  * the -w2 option is the weight for the experimental data score. [default = 0.5]







---

### References:
1. Darty,K., Denise,A. and Ponty,Y. (2009) VARNA: Interactive drawing and editing of the RNA secondary structure. Bioinforma. Oxf. Engl., 25, 1974–197

