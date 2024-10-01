from Functions import run_dna_rna_tools_funcs as seqf

def run_dna_rna_tools(*args):
    # проверка (нет U и T одновременно и иных символов)
    seq_list = list(args[:-1])
    for arg in seq_list:
        t_pos = str(arg).upper().find("T")
        u_pos = str(arg).upper().find("U")
        for i in arg:
            if (i not in "ACGTUacgtu"):
                seq_list.remove(arg)
                break
            if (i in "Uu" and t_pos != -1):
                seq_list.remove(arg)
                break
            if (i == "Tt" and u_pos != -1):
                seq_list.remove(arg)
                break
    # послать в зависимости от последнего аргумента
    if args[-1] == "reverse":
        return seqf.reverse(seq_list)
    if args[-1] == "complement":
        return seqf.complement(seq_list)
    if args[-1] == "reverse_complement":
        return seqf.reverse_complement(seq_list)
    if args[-1] == "transcribe":
        return seqf.transcribe(seq_list)
    if args[-1] == "g_c_percent":
        return seqf.g_c_percent(seq_list)

def filter_fastq(seqs, gc_bound, length_bounds, quality_threshold):
    result = dict()
    return result

