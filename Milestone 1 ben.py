

def dna_count(dna):
    count_A = dna.count('A')
    count_C = dna.count('C')
    count_T = dna.count('T')
    count_G = dna.count('G')
    
    print (count_A,count_C,count_G,count_T)
    return count_A,count_C,count_G,count_T
def dna2rna(dna):

    rna = ''
    for symbol in dna:
        if symbol == 'A':
          rna = rna + 'U'
        elif symbol == 'T':
            rna = rna +'A'
        elif symbol == 'C':
            rna = rna +'G'          
        elif symbol == 'G':
            rna = rna +'C'
  
    print(rna)
    return rna
def reverse_complement(dna):
    reverse_complement = ''
    
    for symbol in dna:
        if symbol == 'A':
            reverse_complement = reverse_complement+'T'
        elif symbol == 'C':
            reverse_complement = reverse_complement+'G'
        elif symbol == 'G':
            reverse_complement = reverse_complement+'C'
        elif symbol == 'T':
            reverse_complement = reverse_complement+'A'
            reverse_dna = dna[::-1]
        
    print(reverse_dna)
    return reverse_dna
dna = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
dna_count(dna)    
dna2rna(dna)
reverse_complement(dna)
