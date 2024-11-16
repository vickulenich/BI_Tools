def get_gc_score(seq: str) -> float:
    '''
    Function gc_score, counts GC bound of sequence

    Args: str

    Returns: float
    '''
    gc_score = (seq.lower().count("g") + seq.lower().count("c"))*100/len(seq)
    return gc_score


def get_quality_score(seq: str) -> float:
    '''
    Function quality_score, counts average phred quality score of sequence

    Args: str

    Returns: float
    '''
    qual_vals = [ord(char)-33 for char in seq]
    return sum(qual_vals)/len(qual_vals)


def read_seq_file(input_fastq: str) -> dict:
    '''
    Function read_sequences, reads initial file and
    forms a dictionary {"sequens_id": ("sequence", "+", "basecall quality")}

    Args: str

    Returns: dict
    '''
    with open(input_fastq, 'r') as input_file:
        seqs = dict()
        lines = input_file.readlines()
        for i in range(0, len(lines), 4):
            seqs[lines[i].rstrip()] = (lines[i + 1].rstrip(),
                                       "+", lines[i + 3].rstrip())
    return seqs
