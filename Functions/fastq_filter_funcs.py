def length_score(seq: str) -> float:
    '''
    Function length_score, counts length of sequence

    Args: str

    Returns: float
    '''
    return len(seq)


def gc_score(seq: str) -> float:
    '''
    Function gc_score, counts GC bound of sequence

    Args: str

    Returns: float
    '''
    gc_score = (seq.lower().count("g") + seq.lower().count("c"))*100/len(seq)
    return gc_score


def quality_score(seq: str) -> float:
    '''
    Function quality_score, counts average phred quality score of sequence

    Args: tuple

    Returns: float
    '''
    qual_vals = list()
    for char in seq:
        qual_vals.append(ord(char)-33)
    return sum(qual_vals)/len(qual_vals)
