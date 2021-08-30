from Sequence_Classes import *
from Utils import *

test = Sequence()
test.generate_rnd_seq(500, "RNA") #this function takes two arguments. sequence length and sequence type. It will generate a random sequence of your chosen length and type (DNA or RNA)

print(test.Seq_Info())

print(test.Nuc_Freq())

print(test.Transcription())

print(test.Reverse_Complement())

print(test.Get_GC_Content())

print(test.Get_Subseq_GC_Content())

print(test.Translation())

print(test.Codon_Use('L'))

for rf in test.Get_Reading_Frames():
    print(rf)

print(test.All_Proteins_From_ORFs())