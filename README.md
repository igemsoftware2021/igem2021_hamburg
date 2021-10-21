# iGEM 2021 Hamburg - LEA

This is LEA (Linker Extraction from Alignments), a command-line tool to extract possible linker sequences for a fusion construct. LEA is a software project that was designed as part of the iGEM 2021 competition. For detailed information about the project and LEA's purpose refer to: https://2021.igem.org/Team:Hamburg/Software

## Installation

### Anaconda

You can use Anaconda to create an environment and install all necessary dependencies. Simply run the following commands:

```bash
conda env create -f environment.yml
conda activate igem-hamburg
```

To deactivate the environment run:
```bash
conda deactivate igem-hamburg
```

Alternatively, you can install the dependencies one by one by running: 
```bash
conda install numpy
conda install requests
conda install -c bioconda mafft
conda install -c bioconda blast
```

### Manual

You can also install the dependencies manually. For NumPy and requests you can use Python's built-in package manager pip by running:

```bash
pip install numpy
pip install requests
```

For installation of MAFFT and BLAST refer to:
- MAFFT: https://mafft.cbrc.jp/alignment/software/
- Blast: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

## Usage

You can execute the full workflow by calling

```bash
python3 lea.py -A <path_to_sequencefile_A> -B <path_to_sequencefile_B> -R <path_to_reference_sequencefile>
```

You can control many different aspects of the workflow with a set of parameters. To get an overview simply run:

```bash
python3 lea.py -h 
```

If you do not want to execute the full workflow you can instead use msa_processor.py to skip the BLAST searches. To Calculate an MSA and extract linker sequences from an existing fasta_file containing the sequence groups run

```bash
python3 msa_processor.py -s <path_to_sequencefile> -m <path_to_msa_outputfile>
```

If you want to use an already existing MSA file for linker extraction you can run

```bash
python3 msa_processor.py -m <path_to_msafile> --skip_msa
```

If you only want to clean an existing MSA file run

```bash
python3 msa_processor.py -m <path_to_msafile> --clean_only
```

## Example

To see an example execution of the workflow run

```bash
sh run_example.sh
```

## License
[MIT](https://choosealicense.com/licenses/mit/)
