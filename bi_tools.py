from Functions import run_dna_rna_tools_funcs as seqf
from Functions import fastq_filter_funcs as fastqf


def run_dna_rna_tools(*args: tuple):
    '''
    Function run_dna_rna_tools, makes transformation
    specified in the last argument

    Args: tuple

    Returns: str or list or float
    '''
    result = list()
    seqs = list(args[:-1])
    for seq in seqs:
        if args[-1] == "reverse":
            result.append(seqf.reverse(seq))
        if args[-1] == "complement":
            result.append(seqf.complement(seq))
        if args[-1] == "reverse_complement":
            result.append(seqf.reverse_complement(seq))
        if args[-1] == "transcribe":
            result.append(seqf.transcribe(seq))
        if args[-1] == "g_c_bound":
            result.append(seqf.get_g_c_bound(seq))
    if len(result) == 1:
        return result[0]
    else:
        return result


def filter_fastq(input_fastq: str, output_fastq: str,
                 gc_bounds: tuple = (0, 100),
                 length_bounds: tuple = (0, 2**32),
                 quality_threshold: int = 0):
    '''
    Function filter_fastq, drops fastq-sequences not meeting specified
    length, gc bound and quality

    Args:
    input_fastq (str)
    output_fastq (str)
    gc_bounds (tuple, default = (0, 100))
    length_bounds (tuple, default = (0, 2**32))
    quality_threshold (int, default = 0)

    Returns: None
    '''
    seqs = fastqf.read_seq_file(input_fastq)
    for key, value in seqs.items():
        if ((gc_bounds[0] <= fastqf.get_gc_score(value[0]) <= gc_bounds[1]) and
                (length_bounds[0] <= len(value[0]) <= length_bounds[1]) and
                (fastqf.get_quality_score(value[2]) >= quality_threshold)):
            with open(output_fastq, 'a') as output_file:
                output_file.write(key + '\n')
                output_file.write('\n'.join(value))
                output_file.write('\n')
