import requests, os, re
from time import sleep
from tempfile import TemporaryFile

def blast_resonse_to_xml_string(rid, rtoe):
    """Takes a Request ID of a Blast search and returns results as a XML-string 

    :param rid: Request ID of the blast search
    :type rid: str
    :param rtoe: Estimated time until response can be expected
    :type rtoe: int
    :raises RuntimeError: If the search returns status FAILED
    :return: XML formatted string of the blast search
    :rtype: str
    """

    sleep(rtoe)

    ready = False
    while not ready:

        check_req = requests.post(
            'https://blast.ncbi.nlm.nih.gov/Blast.cgi',
            data={
                'CMD': 'Get',
                'FORMAT_OBJECT': 'SearchInfo',
                'RID': rid
            }
        )

        # extract status
        status_match = re.search(r'Status=([A-Z]+)', check_req.text)
        if status_match:
          status = status_match.group(1)

        if status == 'READY':
            ready=True
        elif status == 'FAILED':
            raise RuntimeError('Error while executing blast search')

        sleep(5)

    result_req = requests.post(
        'https://blast.ncbi.nlm.nih.gov/Blast.cgi',
        data={
            'CMD': 'Get',
            'FORMAT_TYPE': 'XML',
            'RID': rid
        }
    )

    if result_req.status_code == 200:
        return result_req.content.decode('utf-8')

def response2fasta(response_string, outputfile):
    """saves a blast response file in xml-format into a given fastafile

    :param response_string: Response of a Blast request in xml-format
    :type response_string: str
    :param outputfile: Path to outputfile where the sequences of the blast search are appended
    :type outputfile: str
    """
    response_lines = response_string.splitlines()
    with open(outputfile, 'a') as fastafile:
        for line in response_lines:
            if '<Hit_id>' in line:
                hit_id = re.search(r'(<Hit_id>)(.*)(</Hit_id>)', line).group(2)
                fastafile.write('>' + hit_id)
            if '<Hsp_hseq>' in line:
                hseq = re.search(r'(<Hsp_hseq>)(.*)(</Hsp_hseq>)', line).group(2)
                fastafile.write('\n' + hseq + '\n')

def call_blast_api(sequence, nof_seq):
    """calls the NCBI Blast API and invokes a blastp search with a given sequence

    :param sequence: Protein sequence in FASTA-, textformat or as NCBI accession
    :type sequence: str
    :raises RuntimeError: If the return code of the blast api call is not 200 OK.
    :return: Request ID and estimated time for the blast search by the server
    :rtype: Tuple(str, int)
    """
    invoke_req = requests.post(
        'https://blast.ncbi.nlm.nih.gov/Blast.cgi',
        data={
            'CMD': 'Put',
            'QUERY': sequence,
            'PROGRAM': 'blastp',
            'DATABASE': 'nr',
            'HITLIST_SIZE': nof_seq
        }
    )

    if not invoke_req.status_code == 200:
        raise RuntimeError('Error while invoking blast query.')

    # extract rid and rtoe from response
    rid_match = re.search(r'RID = (([A-Z]|[0-9])+)', invoke_req.text)
    if rid_match:
        rid = rid_match.group(1)

    rtoe_match = re.search(r'RTOE = ([0-9]+)', invoke_req.text)
    if rtoe_match:
        rtoe = int(rtoe_match.group(1))

    return rid, rtoe

def blast_search(sequence_A, sequence_B, output_path, nof_seq):
    """simultaneously conducts a blast search with two sequences and writes the results to a given outputfile

    :param sequence_A: first sequence for blast search in fasta, text format or as NCBI accession
    :type sequence_A: str
    :param sequence_B: second sequence for blast search in fasta, text format or as NCBI accession
    :type sequence_B: str
    :param output_path: Path to outputfile where the sequences of the blast search are appended
    :type output_path: str
    """
    id_A = sequence_A.split("\n")[0]
    id_B = sequence_B.split("\n")[0]
    print(f'Starting web BLAST search with sequence {id_A}')
    rid_A, rtoe_A = call_blast_api(sequence_A, nof_seq)
    print(f'Starting web BLAST search with sequence {id_B}')
    rid_B, rtoe_B = call_blast_api(sequence_B, nof_seq)

    xml_A = blast_resonse_to_xml_string(rid_A, max(rtoe_A, rtoe_B))
    xml_B = blast_resonse_to_xml_string(rid_B, 0)

    response2fasta(xml_A, output_path)
    response2fasta(xml_B, output_path)
    print('Successfully executed web BLAST search. Sequences were added to '+ output_path)

def write_fasta_with_ids(inputfile, out_dir):
    """adds unique ids to a fasta file and removes any unusual characters in 
       and in between the sequences

    :param inputfile: Path to the original fasta file without ids
    :type inputfile: str
    :return: path to the edited fasta file
    :rtype: str
    """
    print('Cleaning reference sequence file for BLAST search.')
    seq_dict = {}
    current_sequence = ''
    current_id = ''

    outputfile = os.path.splitext(os.path.basename(inputfile))[0]
    outputfile = f'{out_dir}/{outputfile}_db.fasta'
    counter = 1
    with open(inputfile, 'r') as db_seq_input, open(outputfile, 'w+') as db_seq_output:
        for line in db_seq_input:
            if '>' in line:
                if current_sequence != '':
                    db_seq_output.write(current_sequence)
                    db_seq_output.write('\n')
                    seq_dict[f'seq{counter}'] = (current_id, current_sequence)
                    current_sequence = ''
                    current_id = ''
                    counter += 1

                current_id += '>seq' + str(counter)
                m = re.findall(r'[A-Za-z0-9()/%.+*"\[\]\-]+', line)
                for word in m:
                    current_id += ' ' + word
                db_seq_output.write(current_id)
                db_seq_output.write('\n')
            else:
                m = re.findall(r'[A-Z]+', line)
                for word in m:
                    current_sequence += word
        db_seq_output.write(current_sequence)

    print('Successfully cleaned reference sequence file. Results were saved to '+outputfile)
    return outputfile, seq_dict

def blastp_local(query_file, subject_file, evalue):
    """executes a blastp search with a given query file in a given subject file

    :param query_file: Path to a fasta file with the query protein sequences
    :type query_file: str
    :param subject_file: Path to a fasta file with sequences that are used as a database for the blast search
    :type subject_file: str
    :param evalue: Maximum evalue for reference sequence hits to be considered for linker extraction
    :type evalue: float
    :return: sequence ids of the hits of the blast search
    :rtype: set()
    """

    # execute blastp with base-sequences as query and known sequences as subject
    blastp_stream = os.popen(
        'blastp -query {} -subject {} -evalue {} -outfmt "7 sseqid evalue"'.format(
                                                            query_file, subject_file, evalue
    ))

    hit_ids = set()

    # extract results, avoid duplicate sequences
    for line in blastp_stream.readlines():
        if '#' not in line:
            seqid = line.split()[0]
            hit_ids.add(seqid)

    return hit_ids

def blast_search_local(queryfile_A, queryfile_B, reference_file, outfile, out_dir, evalue):
    """conducts a blastp search with two query files and a custom subject file 
       through the commandline tool BLAST+ and appends the results to a given outputfile

    :param queryfile_A: Path to query file with sequence(s) in fasta format
    :type queryfile_A: str
    :param queryfile_B: Path to query file with sequence(s) in fasta format
    :type queryfile_B: str
    :param reference_file: Path to subject file with sequences that are searched for hits
    :type reference_file: str
    :param outfile: Path outputfile in fasta format where the hit sequences are appended
    :type outfile: str
    :param evalue: Maximum evalue for reference sequence hits to be considered for linker extraction
    :type evalue: float
    """
    # add unique ids to the given custom fasta sequences
    subject_file, seq_dict = write_fasta_with_ids(reference_file, out_dir)

    print(f'Executing local BLAST searches to filter reference sequences')
    ids_A = blastp_local(queryfile_A, subject_file, evalue)
    ids_B = blastp_local(queryfile_B, subject_file, evalue)
    ids = ids_A.intersection(ids_B)
    write_sequences_from_ids(ids, seq_dict, outfile)
    print('Successfully filtered reference sequences. Remaining sequences were added to '+outfile)

def write_sequences_from_ids(ids, seq_dict, outfile):
    """appends a given outputfile in fasta format with sequences that are given by a list of ids

    :param ids: list of ids
    :type ids: lst(str)
    :param seq_dict: dictionary that translates an id to a sequence
    :type seq_dict: dct[str] = str
    :param outfile: Path outputfile where the sequences are written into in fasta format
    :type outfile: str
    """
    with open(outfile, 'a') as f:
        for seq_id in ids:
            f.write(seq_dict[seq_id][0])
            f.write('\n')
            f.write(seq_dict[seq_id][1])
            f.write('\n')
