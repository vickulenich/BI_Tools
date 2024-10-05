from Functions import run_dna_rna_tools_funcs as seqf
from Functions import fastq_filter_funcs as fastqf


def run_dna_rna_tools(*args: tuple):
    '''
    Function run_dna_rna_tools, makes transformation specified in the last argument

    Args: tuple

    Returns: str or list or float
    '''
    seqs = list(args[:-1])
    for seq in seqs:
        t_position = str(seq).upper().find("T")
        u_position = str(seq).upper().find("U")
        for char in seq:
            if (char not in "ACGTUacgtu"):
                seqs.remove(seq)
                break
            if (char in "Uu" and t_position != -1):
                seqs.remove(seq)
                break
            if (char == "Tt" and u_position != -1):
                seqs.remove(seq)
                break
    if args[-1] == "reverse":
        return seqf.reverse(seqs)
    if args[-1] == "complement":
        return seqf.complement(seqs)
    if args[-1] == "reverse_complement":
        return seqf.reverse_complement(seqs)
    if args[-1] == "transcribe":
        return seqf.transcribe(seqs)
    if args[-1] == "g_c_bound":
        return seqf.g_c_bound(seqs)


def filter_fastq(seqs: dict, gc_bounds: tuple=(0, 100), length_bounds: tuple=(0, 2**32), quality_threshold: int=0):
    '''
    Function filter_fastq, drop fastq-sequences not meeting the specified length, gc bound and quality

    Args: tuple

    Returns: dict
    '''
    for seq in list(seqs):
        if gc_bounds[1] < fastqf.gc_score(seqs[seq][0]) or fastqf.gc_score(seqs[seq][0]) < gc_bounds[0]:
            seqs.pop(seq)
            continue
        elif length_bounds[1] < fastqf.length_score(seqs[seq][0]) or fastqf.length_score(seqs[seq][0]) < length_bounds[0]:
            seqs.pop(seq)
            continue
        elif fastqf.quality_score(seqs[seq][1]) < quality_threshold:
            seqs.pop(seq)
            continue
    return seqs
