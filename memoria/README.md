# Bioinformatics Master degree dissertation

This folder holds my final year project at the [Bioinformatics and Biostatistics](https://estudios.uoc.edu/es/masters-universitarios/bioinformatica-bioestadistica/presentacion) Master degree from the [Universitat Oberta de Catalunya](https://www.uoc.edu).

[PDF](https://github.com/albamgarces/TFM_UOC_AMoya/blob/main/memoria/AlbaMoyaGarces_MemoriaFinal.pdf)(spanish)

Project supervisor: [José Francisco Sánchez Herrero](https://github.com/JFsanchezherrero)

## Abstract
ESKAPE bacteria are thoght to be especially resistant to antibiotics and can be particularly dangerous in people with weakened immune systems. Studies of gene duplication in the genomes of ESKAPE specia and *E. coli* shown differentiated patterns between virulent strains and those that were not, which is a very valuable information that can be useful to struggle cronical diseases. The aim of this TFM is to supplement the previous research works and develops bioinformatics tools in order to carry out a process in a more efficient and intuitive way. Two programs have been implemented on python and R that allow process an annotation file, analize it and carry out all the necessary steps until a annotation file of the duplicated proteins found and plot it. Finally, a small group of selected strains for each ESKAPE species not previously studied has been used as an usage example of the generated tools. Data obtained could be subsequently analyzed in search of relevant patterns or particularities.

### Keywords
microbiology
gene duplication
ESKAPE pathogens

### Experiment Environment:

+ OS: Ubuntu 20.04.1 LTS/x86_64-pc-linux-gnu (64-bit)
+ [Python 3.8.5](https://python.org)
+ [R 4.0.3](https://www.r-project.org/)
+ [RStudio 1.3.1073](https://rstudio.com/)
+ [ncbi-blast+ 2.9.0-2](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

### Preliminary works based

- **Duplicates in *Escherichia coli* genomes**: 

	Gene duplications in the E. coli genome: common themes among pathotypes. Bernabeu M., Sánchez-Herrero JF., Huedo P., Prieto A., Hüttener M., Rozas J. and Juárez A. **BMC Genomics** *2019 20:313*, https://doi.org/10.1186/s12864-019-5683-4

	The source code to replicate these analysis corresponds to the GitHub release of this source code [v1.0](https://github.com/molevol-ub/BacterialDuplicates/releases/tag/v1.0). 

	See an example workflow and additional bioinformatic details and parameters [here](https://github.com/molevol-ub/BacterialDuplicates/blob/master/Ecoli/README.md).

	<br/><br/>

- **Duplicates in Gram positive cocci**:

	Gene Duplications in the Genomes of Staphylococci and Enterococci. Sanchez-Herrero JF., Bernabeu M., Prieto A., Hüttener M. and Juárez A. **Front. Mol. Biosci.** *2020 7:160*. https://doi.org/10.3389/fmolb.2020.00160

	The source code to replicate the analysis corresponds to the GitHub release of this source code [v2.0](https://github.com/molevol-ub/BacterialDuplicates/releases/tag/v2.0). 

	See an example workflow and additional bioinformatic details and parameters [here](https://github.com/molevol-ub/BacterialDuplicates/blob/master/Gram_positive/README.md).

## Usage

### Getting started (Ubuntu)

I am using `pip` and `virtualenv` to create a virtual environment and manage the dependencies.

You will need to have installed `pip` and `virtualenv` to create the virtual environment (more information here: https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/).

- Install pip for python3 (if required):
```
sudo apt install python3-pip
```

- Install python-dev, virtualenv and venv (if required):

```
sudo apt-get install python3-dev
```
```
pip3 install virtualenv
```
```
sudo apt install python3-venv
```

- Create a virtual environment (here called "TFM").
```
python3 -m venv TFM
```

- Activate the environment:
```
source TFM/bin/activate
```
> To deactivate the environment run: `deactivate`. 

### Clone repository

```
git clone https://github.com/albamgarces/TFM_UOC_AMoya.git
cd git/TFM_UOC_AMoya
```

- To install the packages from [`requirements.txt`](https://github.com/albamgarces/memoria/requeriments/pythonenv_requirements.txt):

```
pip3 install -r requirements.txt
```

To get sequences alignment, [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) executables must be downloaded:

```
sudo apt-get install ncbi-blast+
```

To plot the bacterial genomes and duplicate R should be installed and RStudio use are strongly recommended. `Biocircos` library should be installed as shown at [R session](https://github.com/albamgarces/TFM_UOC_AMoya/blob/main/memoria/requeriments/Rsessioninfo) file.

### Duplicate analysis

```
python code/dup_annot.py -h
usage: dupAnnotation [-h] [-a] [-r] [-b] [-bs] [-d] [-e] [-p] [-pi] [-c] [-t]
                     [--pseudo] [-o] [--debug]

Get an annotation file with duplicated protein on genome.

optional arguments:
  -h, --help            show this help message and exit
  -c , --annot_table    Genome annotation .csv file previously analyzed.
  -t , --text_file      Blast raw results text file.
  --pseudo              Wether to use pseudogenes or not
  -o , --out_folder     Results folder
  --debug

Annotation parser options named arguments:
  -a , --annot_file     Annotation file: genbank or GFF.
  -r , --ref_file       Genome references FASTA file.

BLAST options named arguments:
  -b , --blast_folder   BLAST binary folder. **Note blast_folder=/usr/bin is
                        set by default**
  -bs , --bitscore      BLAST bit-score: requires size of a sequence database
                        in which the current match could be found just by
                        chance. **Note bit_score = 50 is set by default**
  -d , --db_name        New database name
  -e , --evalue         BLAST e-value: number of expected hits of similar
                        quality (score) that could be found just by chance.
                        **Note e-value = 1e-05 is set by default**
  -p , --percentage     Percentage of alignment in query. *Note pident = 85 is
                        set by default**
  -pi , --pident        Percentage of similarity in alignment. **Note
                        percentage = 85 is set by default**
```

### Plot results

With RStudio run dup_biocircos.R


## Data

The dataset used in this disseration is shown at the table bellow and available at this [data folder](https://github.com/albamgarces/TFM_UOC_AMoya/tree/main/data) where results are also showed.

![Selected strains to run the tool analysis](https://github.com/albamgarces/TFM_UOC_AMoya/blob/main/memoria/selected_strains.png)

### Data download

Data was downloaded from the [NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/) genome database.

## Output

Files from each strain analysed are also into the correspondent [folder](https://github.com/albamgarces/TFM_UOC_AMoya/tree/main/data).

Genome images are available for *K. pneumoniae* [here](https://github.com/albamgarces/TFM_UOC_AMoya/tree/main/data/Kpneumoniae/biocircos_images)



