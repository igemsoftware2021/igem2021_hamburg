#!/usr/bin/env python3
import numpy as np
import subprocess
import argparse
import json
import os

class MSA_Processor:
    """MSA Processor class. Made to read msa files and extract possible linker sequences.
    """
    def __init__(self,
            sequences_filepath = None,

            msa_filepath = None,
            seq_removal_threshold=5,
            max_gap_width=5,
            gap_perc=0.95,
            num_query_seqs=202,
            out_dir='output_files',
            DEBUG=False):
        """MSA Processor class. Made to read msa files and extract possible linker sequences.

        :param sequences_filepath: Path to a Fasta file. Used as input for MAFFT, defaults to None
        :type sequences_filepath: str, optional
        :param msa_filepath: Path to the MSA file, defaults to None
        :type msa_filepath: str, optional
        :param seq_removal_threshold: Threshold for removal. If the number of non-empty 
                        characters at a given position is below the threshold, the corresponding
                        sequences will be removed, defaults to 1
        :type seq_removal_threshold: int, optional
        :param max_gap_width: Maximal gap witdth that is allowed. Sequences that produce gaps with
                        exceeding width will be removed during cleaning, defaults to 5
        :type max_gap_width: int, optional
        :param gap_perc: Minimum percentage of '-' symbols that a potential linker area needs to
                        show among the query sequences in order to count as a linker area,
                        defaults to 0.95
        :type gap_perc: float, optional
        :param num_query_seqs: The number of query sequences, i.e. all non reference sequences,
                        defaults to 202 (2 query sequences + 100 similar sequences each)
        :type num_query_seqs: int, optional
        :param out_dir: Path to a directory that will hold all output files, defaults to 'output_files'
        :type out_dir: str, optional
        :param DEBUG: Show debug prints, defaults to False
        :type DEBUG: bool, optional
        """
        self.sequences_filepath = sequences_filepath

        self.msa_filepath = msa_filepath
        self.seq_removal_threshold = seq_removal_threshold
        self.max_gap_width = max_gap_width
        self.gap_perc = gap_perc
        self.num_query_seqs = num_query_seqs
        self.DEBUG = DEBUG
        self.identifiers = []
        self.out_dir = out_dir

        self.msa_length = -1
        self.num_seqs = -1
        self.num_ref_seqs = -1

    def write_msa_to_file(self, filename):
        """Write the current msa to a file

        :param filename: Path to a file
        :type filename: str
        """
        with open(filename, 'w') as f:
            i = 0
            for line in self.msa:
                f.write('>seq'+str(i)+'\n')
                f.write(''.join(line))
                f.write('\n')

    def run_msa(self, sequences_filepath = None, output_file = None):
        """Calculate an MSA by using MAFFT

        :param sequences_filepath: Path to the input fasta file, defaults to None
        :type sequences_filepath: str, optional
        :param output_file: Path to the output file, defaults to None
        :type output_file: str, optional
        :raises RuntimeError: If no input file was specified
        :raises RuntimeError: If no output file was specified
        """
        print('Running MAFFT to calculate MSA')
        if sequences_filepath is None:
            sequences_filepath = self.sequences_filepath
        if sequences_filepath is None:
            raise RuntimeError('Error: No fasta input file was specified.')

        if output_file is None:
            output_file = self.msa_filepath
        if output_file is None:
            raise RuntimeError('Error: No MSA file was specified.')

        args = ['mafft',
            '--amino',
            '--quiet',
            '--localpair',
            '--maxiterate 100',
            '--thread 4',
            sequences_filepath,
            '>',
            output_file
        ]

        subprocess.check_call(' '.join(args), shell=True)
        print('Successfully calculated MSA. Saved MSA to '+output_file)

    def load_msa(self, msa_filepath = None):
        """Load an MSA from a file and clean it

        :param msa_filepath: Path to the msa file, defaults to None
        :type msa_filepath: str, optional
        :raises RuntimeError: If no msa file was specified
        """
        if msa_filepath is None:
            msa_filepath = self.msa_filepath
        if msa_filepath is None:
            raise RuntimeError('Error: No MSA file was specified.')

        print('Loading MSA from file '+msa_filepath)

        try:
            self.msa = self._read_msa()
        except BaseException as error:
            print('Error: An error occurred while reading the reading the MSA file.')
            print(str(error))
            return

        print('Cleaning MSA')
        self._update_length_and_num_seqs()
        self.clean_msa()
        self._update_length_and_num_seqs()

        cleaned_path = os.path.splitext(os.path.basename(msa_filepath))[0]
        cleaned_path = f'{self.out_dir}/{cleaned_path}_cleaned.fasta'
        self.write_msa_to_file(cleaned_path)
        print('Successfully loaded and cleaned MSA. Cleaned MSA was saved to '+cleaned_path)

    def clean_msa(self):
        """Data cleanup operations, such as deleting sequences that lead to gaps in the MSA
        """
        if self.DEBUG:
            print(self.msa)
            print()

        change = True
        while change:
            change = False
            remove = self._find_removable_seqs()

            if len(remove) != 0:
                change = True
                self._remove_seqs_and_gaps(remove)
                self._update_length_and_num_seqs()

    def extract_linkers(self):
        """Extract linker sequences from MSA

        :return: Dictionary containing the trimmed query sequences, the linker sequences and their
                corresponding sequence identifiers and the starting position and length of the
                detected linker area in the MSA
        :rtype: dict(str: str)
        """
        print('Extracting linker sequences from MSA')
        if self.msa is None:
            raise RuntimeError('Error: No MSA was loaded.')

        linker_start, linker_length = self._extract_linker_region()

        if linker_start == -1:
            raise RuntimeError('Error: No linker area could be detected')

        linker_msa = self.msa[self.num_query_seqs:, linker_start:linker_start+linker_length]

        new_sequence_A = ''.join(self.msa[0, 0:linker_start])
        new_sequence_A = new_sequence_A.replace('-','')

        new_sequence_B = ''.join(self.msa[1, linker_start+linker_length:])
        new_sequence_B = new_sequence_B.replace('-','')

        linkers = [''.join(seq) for seq in linker_msa]
        linkers = [linker.replace('-', '') for linker in linkers]

        linker_dict = {self.identifiers[i+self.num_query_seqs]: l for i, l in enumerate(linkers)}
        result = {'trimmed_sequence_A': new_sequence_A,
                'trimmed_sequence_B': new_sequence_B,
                'linker_start_in_msa': linker_start,
                'linker_length_in_msa': linker_length,
                'linker_sequences': linker_dict}
        print('Successfully extracted linker sequences. Result: \n'+json.dumps(result, indent = 4))
        return result

    def _read_msa(self):
        """Read MSA from file

        :return: MSA
        :rtype: numpy.ndarray, shape: (x,x)
        """
        with open(self.msa_filepath, 'r') as msa_file:
            msa_lines = []
            new_seq = False
            for line in msa_file.readlines():
                if line.startswith('>'):
                    self.identifiers.append(line[1:].rstrip())
                    new_seq = True
                    continue
                if new_seq:
                    msa_lines.append(line.rstrip())
                    new_seq = False
                else:
                    msa_lines[-1] = msa_lines[-1] + line.rstrip()

        if self.DEBUG:
            print(self.identifiers)

        msa = [list(l[:-1]) for l in msa_lines]
        return np.array(msa)

    def _update_length_and_num_seqs(self):
        """Update the msa length and number of sequence variables after removal of sequences
        """
        self.msa_length = len(self.msa[0])
        self.num_seqs = len(self.msa)
        self.num_ref_seqs = self.num_seqs - self.num_query_seqs

    def _find_removable_seqs(self):
        """Identify the sequences that should be removed as they lead to gaps

        :return: List of sequence indices that should be removed
        :rtype: list
        """
        remove = set()
        suspicious = set()
        for pos in range(self.msa_length-self.max_gap_width):
            msa_window = self.msa[:,pos:pos+self.max_gap_width]
            for seq_idx in range(self.num_seqs):
                if not np.all(msa_window[seq_idx] == '-'):
                    suspicious.add(seq_idx)
            if len(suspicious) <= self.seq_removal_threshold:
                remove = remove.union(suspicious)
            suspicious = set()
        return list(remove)

    def _remove_seqs_and_gaps(self, remove):
        """Remove sequences and resulting gaps from MSA

        :param remove: List of sequence indices that should be removed
        :type remove: list
        """
        msa = np.delete(self.msa, remove, 0)
        self.identifiers = np.delete(self.identifiers, remove)
        self.num_query_seqs -= np.count_nonzero(np.array(remove) < self.num_query_seqs)

        if self.DEBUG:
            print(msa)
            print()

        gaps = np.nonzero([np.all(msa[:,pos]=='-') for pos in range(self.msa_length)])
        msa = np.delete(msa, gaps, 1)

        if self.DEBUG:
            print(msa)
            print()
            print()
            print()

        self.msa = msa

    def _extract_linker_region(self):
        """Identify beginning and length of the most probable linker region

        :return: Tuple containing the linker starting position and length
        :rtype: Tuple(int,int)
        """
        possible_linkerpositions = [
            # For each column, there must be more than <gap_perc> percent of '-' symbols in the query sequence groups
            np.count_nonzero(self.msa[:self.num_query_seqs,pos]=='-') > self.num_query_seqs*self.gap_perc
            for pos in range(self.msa_length)
        ]

        linker_start = -1
        longest_linker = 0
        current_linker = 0
        last_was_linker = False
        for i, is_linker in enumerate(possible_linkerpositions):
            if last_was_linker and is_linker:
                current_linker += 1
            elif last_was_linker and not is_linker:
                if longest_linker < current_linker:
                    longest_linker = current_linker
                    linker_start = i - current_linker
                current_linker = 0
                last_was_linker = False
            elif not last_was_linker and is_linker:
                current_linker += 1
                last_was_linker = True

        return linker_start, longest_linker

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate MSAs and extract possible linker sequences')
    parser.add_argument('-s', '--sequences_filepath', required=False, action='store',
            help="Path to a Fasta file containing sequences for MSA calculation")
    parser.add_argument('-m', '--msa_filepath', required=False, action='store',
            help="Path to where the MSA file will be stored or to an existing MSA file. If you provide "+\
            "an existing MSA make sure the first two sequences are the query sequences A and B and that "+\
            "you check the --skip_msa option")
    parser.add_argument('--skip_msa', action='store_true',
            help="Set to True if no new msa should be calculated")
    parser.add_argument('--clean_only', action='store_true',
            help="Set to True if you only want to clean an MSA file")
    parser.add_argument('-o', '--out_dir', action='store', required=False, default='output_files',
            help="Path to a directory that will hold all output files")
    parser.add_argument('-srt', '--seq_removal_threshold', default=5, action='store', type=int,
            help="Threshold for removal of sequences during MSA cleanup. If the number of non-empty "+\
            "characters at a given position is below the threshold, the corresponding "+\
            "sequences will be removed")
    parser.add_argument('-mgw', '--max_gap_width', default=5, action='store', type=int,
            help="Maximal gap witdth that is allowed. Sequences that produce gaps with "+\
            "exceeding width will be removed during cleaning")
    parser.add_argument('-gp', '--gap_perc', default=0.95, action='store', type=float,
            help="Minimum percentage of '-' symbols that a potential linker area needs to "+\
            "show among the query sequences in order to count as a linker area")
    parser.add_argument('-n', '--num_query_seqs', default=202, action='store', type=int,
            help="The number of query sequences, i.e. all non reference sequences, "+\
            "defaults to 202 (2 query sequences + 100 similar sequences each)")
    args = vars(parser.parse_args())

    if args['sequences_filepath'] is not None and not os.path.exists(args['sequences_filepath']):
        raise RuntimeError('File does not exist: ' + args['sequences_filepath'])

    if args['msa_filepath'] is not None and not os.path.exists(args['msa_filepath']):
        raise RuntimeError('File does not exist: ' + args['msa_filepath'])

    if args['out_dir'] is not None and not os.path.exists(args['out_dir']):
        os.makedirs(args['out_dir'])

    skip_msa = args.pop('skip_msa')
    clean_only = args.pop('clean_only')

    processor = MSA_Processor(**args)

    try:
        if not skip_msa and not clean_only:
            processor.run_msa()
        processor.load_msa()
        if clean_only:
            exit(0)
        result = processor.extract_linkers()
    except RuntimeError as error:
        print(str(error))

    res_file = f'{args["out_dir"]}/linker_sequences.txt'
    with open(res_file, 'w') as f:
        f.write(json.dumps(result, indent = 4))

    print('Done. Results were saved to '+res_file)
