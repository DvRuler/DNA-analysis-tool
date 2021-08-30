from Nucleotides_Codons import *
from collections import Counter
import random

class Sequence:
    #Initialise and validate sequence
    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__Seq_Validate()
        assert self.is_valid, f"Given sequence is not a valid {self.seq_type} sequence"

    def __Seq_Validate(self):
        #Check the sequence to make sure it is a valid DNA/RNA string
        return set(Nucleotides[self.seq_type]).issuperset(self.seq)

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        #Generate a random DNA sequence, provided the length
        seq = ''.join([random.choice(Nucleotides[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def Seq_Info(self):
        #Returns sequence information
        return f'[Label]: {self.label}\n[Sequence]: {self.seq}\n[Type]: {self.seq_type}\n[Length]: {len(self.seq)}'

    def Seq_Type(self):
        #Returns sequence type: DNA/RNA/"""
        return self.seq_type

    def Nuc_Freq(self):
        #Count nucleotide frequency within DNA seq
        return dict(Counter(self.seq))

    def Transcription(self):
        #Returns RNA seq from given DNA sequence
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a valid DNA sequence"    

    def Reverse_Complement(self):
        #Returns reverse complement of given sequence
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def Get_GC_Content(self):
        return round((self.seq.count("C") + self.seq.count("G"))/len(self.seq) * 100)

    def Get_Subseq_GC_Content(self, k=5):
        GC_res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            GC_res.append(round((subseq.count("C") + subseq.count("G"))/len(subseq) * 100))
        return GC_res

    def Translation(self, start_pos=0):
        #Translates DNA or RNA seq to amino acid seq
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(start_pos, len(self.seq) -2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(start_pos, len(self.seq) -2, 3)]

    def Codon_Use(self, amino_acid):
        #Returns frequency of each amino acid encoding codon used in sequence
        templist = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == amino_acid:
                    templist.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == amino_acid:
                    templist.append(self.seq[i:i + 3])

        freqDict = dict(Counter(templist))
        total = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / total, 2)
        return freqDict      

    def Get_Reading_Frames(self):
        #Returns open reading frames from a sequence as well as its reverse complement
        frames = []
        frames.append(self.Translation(0))
        frames.append(self.Translation(1))
        frames.append(self.Translation(2))
        tmp_Seq = Sequence(self.Reverse_Complement(), self.seq_type)
        frames.append(tmp_Seq.Translation(0))
        frames.append(tmp_Seq.Translation(1))
        frames.append(tmp_Seq.Translation(2))
        del tmp_Seq
        return frames

    def ORF_Proteins(self, AA_seq):
        #Returns all possible proteins from a sequence
        current_prot = []
        proteins = []
        for AA in AA_seq:
            if AA == "_":
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                if AA == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += AA     
        return proteins

    def All_Proteins_From_ORFs(self, StartReadPos = 0, EndReadPos = 0, ordered = False):
        #Returns all possible proteins from generated reading frames
        if EndReadPos > StartReadPos:
            tmp_Seq = Sequence(self.seq[StartReadPos: EndReadPos], self.seq_type)
            rfs = tmp_Seq.Get_Reading_Frames()
        else:
            rfs = self.Get_Reading_Frames()
            
        res = []
        for rf in rfs:
            prots = self.ORF_Proteins(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key = len, reverse = True)
        return res