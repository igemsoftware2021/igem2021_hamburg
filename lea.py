#!/usr/bin/env python3
import os
from tempfile import TemporaryDirectory
from msa_processor import MSA_Processor
import argparse
import json

from blast import blast_search, blast_search_local

def extract_sequences(path_A, path_B, msa_input_file):
    """Extract sequences from input files and write them to the MSA input file for later use

    :param path_A: Path to fasta file of sequence A
    :type path_A: str
    :param path_B: Path to fasta file of sequence B
    :type path_B: str
    :param msa_input_file: Path to the MSA input file used later
    :type msa_input_file: str
    :return: Sequence A and Sequence B
    :rtype: Tuple(str, str)
    """
    with open(path_A, 'r') as f_A,\
            open(path_B, 'r') as f_B,\
            open(msa_input_file, 'w+') as f_MSA:
        sequence_A = f_A.read()
        sequence_B = f_B.read()
        f_MSA.write(sequence_A)
        f_MSA.write(sequence_B)
    return sequence_A, sequence_B

def main(sequence_A_path, sequence_B_path, ref_sequences_path, out_dir='output_files', evalue=0.01,
        seq_removal_threshold=5, max_gap_width=5, gap_perc=0.95, num_seqs_per_group=100):
    """Main workflow. Build sequence groups with blast and extract possible linker sequences from
    reference sequences. The goal is to find Linker sequences L that can connect sequences A and B
    to an operational fusion protein of the form A-L-B

    :param sequence_A_path: Path to fasta file of sequence A
    :type sequence_A_path: str
    :param sequence_B_path: Path to fasta file of sequence B
    :type sequence_B_path: str
    :param ref_sequences_path: Path to fasta file with reference sequences
    :type ref_sequences_path: str
    :param out_dir: Path to a directory that will hold all output files, defaults to 'output_files'
    :type out_dir: str, optional
    :param evalue: Maximal evalue for the local blast search. Only reference sequences that match
                    both query sequences with a lower evalue will be considered in the linker
                    extraction
    :type evalue: float
    :param seq_removal_threshold: Threshold for removal. If the number of non-empty 
                    characters at a given position is below the threshold, the corresponding
                    sequences will be removed, defaults to 5
    :type seq_removal_threshold: int, optional
    :param max_gap_width: Maximal gap witdth that is allowed. Sequences that produce gaps with
                    exceeding width will be removed during cleaning, defaults to 5
    :type max_gap_width: int, optional
    :param gap_perc: Minimum percentage of '-' symbols that a potential linker area needs to
                    show among the query sequences in order to count as a linker area,
                    defaults to 0.95
    :type gap_perc: float, optional
    :param num_query_seqs: The number of sequences to take from the blast search of the query 
                    sequences to form a group, defaults to 100
    :type num_query_seqs: int, optional
    """
    basename_A = os.path.splitext(os.path.basename(sequence_A_path))[0]
    basename_B = os.path.splitext(os.path.basename(sequence_B_path))[0]
    basename_R = os.path.splitext(os.path.basename(ref_sequences_path))[0]

    name_ext = f'{basename_A}_{basename_B}_{basename_R}'

    msa_input_path =\
        f'{out_dir}/combined_sequences_{name_ext}.fasta'

    sequence_A, sequence_B = extract_sequences(sequence_A_path,
                                                sequence_B_path,
                                                msa_input_path)

    blast_search(sequence_A, sequence_B, msa_input_path, num_seqs_per_group)

    blast_search_local(sequence_A_path, sequence_B_path, ref_sequences_path, msa_input_path,
                    out_dir, evalue)

    msa_processor = MSA_Processor(sequences_filepath=msa_input_path, 
                                msa_filepath=f'{out_dir}/msa_{name_ext}.fasta',
                                seq_removal_threshold=seq_removal_threshold,
                                max_gap_width=max_gap_width,
                                gap_perc=gap_perc,
                                num_query_seqs=2*num_seqs_per_group+2)

    try:
        msa_processor.run_msa()
        msa_processor.load_msa()
        result = msa_processor.extract_linkers()
    except RuntimeError as error:
        print(str(error))

    result = json.dumps(result, indent = 4)

    res_file = f'{out_dir}/linker_sequences_{name_ext}.txt'
    with open(res_file, 'w') as f:
        f.write(result)

    print('Done. Results were saved to '+res_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract linker sequences from reference proteins. '+\
        'The goal is to find Linker sequences L that can connect sequences A and B '+\
        'to create an operational fusion protein of the form A-L-B')
    parser.add_argument('-A', '--sequence_A_path', required=True, action='store',
            help="Path to a Fasta file containing the first protein's amino acid sequence")
    parser.add_argument('-B','--sequence_B_path', required=True, action='store',
            help="Path to a Fasta file containing the second protein's amino acid sequence")
    parser.add_argument('-R', '--ref_sequences_path', required=True, action='store',
            help="Path to a Fasta file containing amino acid sequences of refecence proteins")
    parser.add_argument('-o', '--out_dir', default='output_files', action='store',
            help="Path to a directory that will hold all output files")
    parser.add_argument('-e', '--evalue', default=0.01, action='store',
            help="Maximal evalue for the local blast search. Only reference sequences that match "+\
            "both query sequences with a lower evalue will be considered in the linker extraction")
    parser.add_argument('-srt', '--seq_removal_threshold', default=5, action='store',
            help="Threshold for removal of sequences during MSA cleanup. If the number of non-empty "+\
            "characters at a given position is below the threshold, the corresponding "+\
            "sequences will be removed")
    parser.add_argument('-mgw', '--max_gap_width', default=5, action='store', type=int,
            help="Maximal gap witdth that is allowed. Sequences that produce gaps with "+\
            "exceeding width will be removed during cleaning")
    parser.add_argument('-gp', '--gap_perc', default=0.95, action='store',
            help="Minimum percentage of '-' symbols that a potential linker area needs to "+\
            "show among the query sequences in order to count as a linker area")
    parser.add_argument('-n', '--num_seqs_per_group', default=100, action='store',
            help="The number of sequences to take from the blast search of the query sequences to form a group")
    args = vars(parser.parse_args())

    if not os.path.exists(args['sequence_A_path']):
        raise RuntimeError('File does not exist: ' + args['sequence_A_path'])

    if not os.path.exists(args['sequence_B_path']):
        raise RuntimeError('File does not exist: ' + args['sequence_B_path'])

    if not os.path.exists(args['ref_sequences_path']):
        raise RuntimeError('File does not exist: ' + args['ref_sequences_path'])

    if 'out_dir' in args.keys() and not os.path.exists(args['out_dir']):
        os.makedirs(args['out_dir'])

    main(**args)
