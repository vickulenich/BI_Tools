def reverse(*args):  # разворачивает
    result = list()
    for arg in list(*args):
        rev_seq = str(arg[::-1])
        result.append(rev_seq)
    if len(result) == 1:
        return result[0]
    else:
        return result
    
def complement(*args):  # заменяет на комплементарную
    result = list()
    init_line = "ACGTUacgtu"
    exit_line = "TGCAAtgcaa"
    for arg in list(*args):
        tab = str(arg).maketrans(init_line, exit_line)
        compl_seq = str(arg).translate(tab)
        result.append(compl_seq)
    if len(result) == 1:
        return result[0]
    else:
        return result
    
def reverse_complement(*args):  # разворачивает и заменяет на комплементарную
    result = list()
    for arg in list(args):
        rev_seq = reverse(arg)
        rev_compl_seq = complement(rev_seq)
        result.append(''.join(rev_compl_seq))
    if len(result) == 1:
        return result[0]
    else:
        return result
    
def transcribe(*args):  # ДНК - вернет РНК, РНК - саму себя
    result = list()
    for arg in list(*args):
        if str(arg).upper().find("U") != -1:
            result.append(arg)
        else:
            trans_seq = arg.replace("T", "U").replace("t", "u")
            result.append(trans_seq)
    if len(result) == 1:
        return result[0]
    else:
        return result
    
def g_c_percent(*args):  # считает ГЦ-состав последовательности
    result = list()
    for arg in list(*args):
        percent = (str(arg).lower().count("g") +
                   str(arg).lower().count("c"))*100/len(arg)
        result.append(round(percent, 2))
    if len(result) == 1:
        result = int(result)
    else:
        return result
