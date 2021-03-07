
def dna_count(dna):
    dna = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
    dna = dna.upper()
    count_A = dna.count('A')
    count_C = dna.count('C')
    count_T = dna.count('T')
    count_G = dna.count('G')
    return count_A,count_C,count_G,count_T
    print (count_A,count_C,count_G,count_T)
def dna2rna(dna):
    dna ='TGCA'
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
    return rna
    print(rna)
def reverse_complement(dna):
    reverse_complement = ''
    dna ='TGCA'
    for symbol in dna:
        if symbol == 'A':
            reverse_complement = reverse_complement+'T'
        elif symbol == 'C':
            reverse_complement = reverse_complement+'G'
        elif symbol == 'G':
            reverse_complement = reverse_complement+'C'
        elif symbol == 'T':
            reverse_complement = reverse_complement+'A'
    return(reverse_complement)
    print(reverse_complement)