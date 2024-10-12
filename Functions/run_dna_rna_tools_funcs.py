def reverse(seq: str) -> str:
    '''
    Function reverse, transforms initial sequence to reversed

    Args: str

    Returns: str
    '''
    return str(seq[::-1])


def complement(seq: str) -> str:
    '''
    Function complement, transforms initial sequence to complemented

    Args: str

    Returns: str
    '''
    init_line = "ACGTUacgtu"
    exit_line = "TGCAAtgcaa"
    tab = str(seq).maketrans(init_line, exit_line)
    return str(seq).translate(tab)


def reverse_complement(seq: str) -> str:
    '''
    Function reverse_complement, transforms
    initial sequence to reversed and complemented

    Args: str

    Returns: str
    '''
    rev_seq = reverse(seq)
    rev_compl_seq = complement(rev_seq)
    return ''.join(rev_compl_seq)


def transcribe(seq: str) -> str:
    '''
    Function transcribe, transforms initial sequence
    to RNA if DNA or returns inintial sequence if RNA

    Args: str

    Returns: str
    '''
    if str(seq).upper().find("U") != -1:
        return seq
    else:
        return seq.replace("T", "U").replace("t", "u")


def get_g_c_bound(seq: str) -> float:
    '''
    Function g_c_bound, counts GC bound of sequence

    Args: str

    Returns: float
    '''
    gc_bound = (str(seq).lower().count("g")
                + str(seq).lower().count("c"))*100/len(seq)
    return gc_bound
