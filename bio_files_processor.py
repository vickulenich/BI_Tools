def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str = None):
    '''
    Function convert_multiline_fasta_to_oneline,

    Args:
    input_fasta (str)
    output_fasta (str, optional)

    Returns: None
    '''
    with open(input_fasta, 'r') as input_file:
        fasta_data = input_file.read()
    oneline_seqs = list()
    for seq in fasta_data.split('>')[1:]:
        lines = seq.split('\n')
        oneline_seqs.append('>' + lines[0] + '\n' +
                            ''.join(lines[1:]).replace('\n', '') + '\n')
    oneline_fasta = ''.join(oneline_seqs)
    if output_fasta:
        with open(output_fasta, 'w') as output_file:
            output_file.write(oneline_fasta)
    else:
        print(oneline_fasta)


def parse_blast_output(input_file: str, output_file: str):
    '''
    Function convert_multiline_fasta_to_oneline,

    Args:
    input_file (str)
    output_file (str)

    Returns: None
    '''
#    with open(input_file, 'r') as input_file:
#        blast_data = input_file.read()
#
#    gene_sequences = list()
#
#    with open(output_file, 'w') as output_file:
#        for protein in proteins:
#            output_file.write(protein + '\n')


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list,
                                   output_fasta: str, n_before: int = 1,
                                   n_after: int = 1):
    '''
    Function convert_multiline_fasta_to_oneline,

    Args:
    input_gbk (str)
    genes (list)
    output_fasta (str)
    n_before (int, default = 1)
    n_after (int, default = 1)

    Returns: None
    '''
#    with open(input_gbk, 'r') as input_file:
#        gbk = input_file.read()
#
#    gene_seqs = list()
#
#    with open(output_fasta, 'w') as output_file:
#        output_file.write(''.join(gene_seqs))
