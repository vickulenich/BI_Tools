from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import GC


class BiologicalSequence(ABC):
    def __init__(self, sequence):
        self.sequence = sequence

    @abstractmethod
    def get_length(self) -> int:
        pass

    @abstractmethod
    def get_subsequence(self, start: int, end: int) -> str:
        pass

    @abstractmethod
    def print_sequence(self) -> None:
        pass

    @abstractmethod
    def check_alphabet(self, sequence: str) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)

    def get_length(self) -> int:
        return len(self.sequence)

    def get_subsequence(self, start, end) -> str:
        return self.sequence[start:end+1]

    def print_sequence(self) -> None:
        print(self.sequence)

    def complement(self) -> str:
        '''
        Function complement, transforms initial sequence to complemented

        Args: str

        Returns: str
        '''
        init_line = "ACGTUacgtu"
        exit_line = "TGCAAtgcaa"
        tab = str(self.sequence).maketrans(init_line, exit_line)
        return str(self.sequence).translate(tab)

    def reverse(self) -> str:
        '''
        Function reverse, transforms initial sequence to reversed

        Args: str

        Returns: str
        '''
        return str(self.sequence[::-1])

    def reverse_complement(self) -> str:
        '''
        Function reverse_complement, transforms initial sequence to reversed and complemented

        Args: str

        Returns: str
        '''
        self.sequence = self.reverse()
        rev_compl_seq = self.complement()
        return ''.join(rev_compl_seq)

    def get_g_c_score(self) -> float:
        '''
        Function g_c_bound, counts GC bound of sequence

        Args: str

        Returns: float
        '''
        gc_score = (str(self.sequence).lower().count("g")
                    + str(self.sequence).lower().count("c"))*100/self.get_length()
        return gc_score

    def check_alphabet(self) -> bool:
        valid_nucleotides = set("ACGTUacgtu")
        return all(nucleotide in valid_nucleotides for nucleotide in self.sequence)


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)

    def transcribe(self) -> str:
        '''
        Function transcribe, transforms initial DNA sequence to RNA

        Args: str

        Returns: str
        '''
        return self.sequence.replace("T", "U").replace("t", "u")

    def check_alphabet(self) -> bool:
        valid_nucleotides = set("ACGTacgt")
        return all(nucleotide in valid_nucleotides for nucleotide in self.sequence)


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)

    def check_alphabet(self) -> bool:
        valid_nucleotides = set("ACGUacgu")
        return all(nucleotide in valid_nucleotides for nucleotide in self.sequence)


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)

    def get_length(self) -> int:
        return len(self.sequence)

    def get_subsequence(self, start: int, end: int) -> str:
        return self.sequence[start:end+1]

    def print_sequence(self) -> None:
        print(self.sequence)

    def get_aa_percentage(self) -> None:
        polar_count = 0
        nonpolar_count = 0

        polar_acids = set("GAVLIPFW")
        nonpolar_acids = set("STCMNQYDEKRN")

        polar_count = sum(1 for aa in self.sequence if aa in polar_acids)
        nonpolar_count = sum(1 for aa in self.sequence if aa in nonpolar_acids)

        print(f"The amino acid sequence has {polar_count/self.get_length()*100}% of polar acids and {nonpolar_count/self.get_length()*100}% of non-polar acids.")

    def check_alphabet(self) -> bool:
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        return all(amino_acid in valid_amino_acids for amino_acid in self.sequence)


def filter_fastq(input_fastq: str, output_fastq: str,
                 gc_bounds: tuple = (0, 100),
                 length_bounds: tuple = (0, 2**32),
                 quality_threshold: int = 0):
    '''
    Function filter_fastq, drops fastq-sequences not meeting specified length, gc bound and quality

    Args:
    input_fastq (str): path to the input FASTQ file.
    output_fastq (str): path to the output FASTQ file.
    gc_bounds (tuple, default = (0, 100)): lower and upper bounds for GC content.
    length_bounds (tuple, default = (0, 2**32)): lower and upper bounds for sequence length.
    quality_threshold (int, default = 0): min quality score.

    Returns: None
    '''
    for record in SeqIO.parse(input_fastq, "fastq"):
        seq = str(record.seq)
        gc_content = GC(seq)
        quality_score = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
        if ((gc_bounds[0] <= gc_content <= gc_bounds[1]) and
                (length_bounds[0] <= len(seq) <= length_bounds[1]) and
                (quality_score >= quality_threshold)):
            with open(output_fastq, 'a') as output_file:
                SeqIO.write(record, output_file, "fastq")
