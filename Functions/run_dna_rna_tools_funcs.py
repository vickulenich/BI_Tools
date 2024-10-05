def reverse(*args: tuple):
    '''
    Function reverse, transforms initial sequences to reversed

    Args: tuple

    Returns: str or list
    '''
    rev_seqs = list()
    for seq in list(*args):
        rev_seq = str(seq[::-1])
        rev_seqs.append(rev_seq)
    if len(rev_seqs) == 1:
        return rev_seqs[0]
    else:
        return rev_seqs


def complement(*args: tuple):
    '''
    Function complement, transforms initial sequences to complemented

    Args: tuple

    Returns: str or list
    '''
    compl_seqs = list()
    init_line = "ACGTUacgtu"
    exit_line = "TGCAAtgcaa"
    for seq in list(*args):
        tab = str(seq).maketrans(init_line, exit_line)
        compl_seq = str(seq).translate(tab)
        compl_seqs.append(compl_seq)
    if len(compl_seqs) == 1:
        return compl_seqs[0]
    else:
        return compl_seqs


def reverse_complement(*args: tuple):
    '''
    Function reverse_complement, transforms each of initial sequences to reversed and complemented

    Args: tuple

    Returns: str or list
    '''
    rev_compl_seqs = list()
    for seq in list(args):
        rev_seq = reverse(seq)
        rev_compl_seq = complement(rev_seq)
        rev_compl_seqs.append(''.join(rev_compl_seq))
    if len(rev_compl_seqs) == 1:
        return rev_compl_seqs[0]
    else:
        return rev_compl_seqs


def transcribe(*args: tuple):
    '''
    Function transcribe, transforms each of initial sequences to RNA if DNA or returns inintial sequence if RNA

    Args: tuple

    Returns: str or list
    '''
    transcr_seqs = list()
    for seq in list(*args):
        if str(seq).upper().find("U") != -1:
            transcr_seqs.append(seq)
        else:
            trans_seq = seq.replace("T", "U").replace("t", "u")
            transcr_seqs.append(trans_seq)
    if len(transcr_seqs) == 1:
        return transcr_seqs[0]
    else:
        return transcr_seqs


def g_c_bound(*args: tuple):
    '''
    Function g_c_bound, counts GC bound of each sequence

    Args: tuple

    Returns: int or list
    '''
    gc_bounds = list()
    for seq in list(*args):
        bound = (str(seq).lower().count("g") +
                 str(seq).lower().count("c"))*100/len(seq)
        gc_bounds.append(round(bound, 2))
    if len(gc_bounds) == 1:
        gc_bounds = int(gc_bounds)
    else:
        return gc_bounds
