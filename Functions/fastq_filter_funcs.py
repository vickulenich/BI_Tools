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

    Args: tuple

    Returns: float
    '''
    qual_vals = [ord(char)-33 for char in seq]
    return sum(qual_vals)/len(qual_vals)
