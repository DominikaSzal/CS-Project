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

  a = dna_list[0]
  b = dna_list[1]
  c = dna_list[2]

  GC_count0 = a.count('C') + a.count('G')
  GC_count1 = b.count('C') + b.count('G')
  GC_count2 = c.count('C') + c.count('G')

  if GC_count0 > GC_count1 and GC_count0 > GC_count2:
      return (0, (GC_count0/len(dna_list[0]))*100)
  elif GC_count1 > GC_count0 and GC_count1 > GC_count2:
      return (1, (GC_count1/len(dna_list[1]))*100)
  elif GC_count2 > GC_count0 and GC_count2 > GC_count1:
      return (2, (GC_count2/len(dna_list[2]))*100)

def locate_substring(dna_snippet, dna): 
        indexes = [i for i in range(len(dna_snippet)) if dna_snippet.startswith(dna, i)] 
        return indexes

print(locate_substring("GATATATGCATATACTT","ATAT"))

def hamming_dist(dna1, dna2):
    count = 0
    d = 0
    while(d<len(dna1)):
        if dna1[d] != dna2[d]:
            count += 1
        d += 1
    return count

def rna2codon(triplet):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y',                                'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CGU': 'R',        'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    allowed_codons = set('ACGU')
    tripletUP = triplet.upper()
    for character in tripletUP:
        if character not in allowed_codons:
            return('')
    return_val = genetic_code.get(tripletUP, '')

    return return_val

def rna2codons(rna):
    amino_string = ''
    for i in range(0,int(len(rna)/3)):
        c = rna[3*i:3*i+3]
        amino = rna2codon(c)
        amino_string = amino_string + amino
    return amino_string

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
