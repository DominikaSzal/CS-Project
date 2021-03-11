def dna_count(dna):
    count_A = dna.count('A')
    count_C = dna.count('C')
    count_T = dna.count('T')
    count_G = dna.count('G')

    d = {'A':count_A, 'C':count_C, 'T':count_T, 'G':count_G}
    print (d)
    return d

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

   for symbol in dna:
       if symbol == 'A':
           dna = dna + 'T'
       elif symbol == 'T':
           dna = dna + 'A'
       elif symbol == 'C':
           dna = dna + 'G'
       elif symbol == 'G':
           dna = dna + 'C'
   reverse_dna = dna[::-1]

   return reverse_dna

def mendels_law(hom,het,rec):
    population = hom + het + rec
    pop2 = population - 1
    
    het_pop = het/population
    rec_pop = rec/population

    
    het_pop2 = (het - 1)/pop2
    het_pop2b = het/pop2
    rec_pop3 = (rec - 1)/pop2
    rec_pop3b = rec/pop2

    return 1 - rec_pop*rec_pop3 - rec_pop*het_pop2b*0.5 - het_pop*rec_pop3b*0.5 - het_pop*het_pop2*0.25



def fibonacci_rabbits(n,k):
    F1,F2 = 1,1
    for i in range(n-1):
        F2,F1 = F1,F1+(F2*k)
    
    return F2


dna_list = ['CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG','CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC','CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT']

def GC_content(dna_list):
    percentages = []
    for index in range(len(dna_list)):
        C = dna_list[index].count('C')
        G = dna_list[index].count('G')
        total = len(dna_list[index])
        count = C + G
        GC_percent = (count/total)*100
        percentages.append(CG_percent)
    highest_percentage = max(percentages)
    index_highest_percentage = percentages.index(highest_percentage)
    return ((index_highest_percentage,highest_percentage))

def locate_substring(dna_snippet, dna): 
        indexes = [i for i in range(len(dna_snippet)) if dna_snippet.startswith(dna, i)] 
        return indexes



def hamming_dist(dna1, dna2):
    count = 0
    d = 0
    while(d<len(dna1)):
        if dna1[d] != dna2[d]:
            count += 1
        d += 1
    return count

 def rna2codon_i(rnas):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    amino = ''
    if rnas in genetic_code:
        amino - genetic_code[str(rna)]
    else:
        amino = 'Invalid'
    return amino
    

def rna2codon(rna):
    amino = ''
    for i in range(0,int(len(rna)/3)):
        if rna2codon_i(rna[3*i:3*i_3]) == '*':
            return amino
        else:
            amino = amino + rna2codon_i(rna[3*i:3*i+3])
    return amino

def count_dom_phenotype(genotypes):
    offspring = (genotypes[0] + genotypes[1] + genotypes[2])*2 + genotypes[3]*1.5 + genotypes[4] 
    return offspring

def source_rna(protein):
  count = 1
  for key in genetic_code.keys():
      if genetic_code[key] in protein:
          count += 1
      elif genetic_code[key]=='*':
          count += len(protein)
  return count % 100000

def splice_rna(dna, intron_list):
    for x in intron_list:
        for i in dna:
            dna = dna.replace(x,"")
    
    rna = dna2rna(dna)

    #see campuswire post #867
    for item in rna:
            rna = rna.replace("T","U")

    protein = rna2codon(rna)

    return protein
