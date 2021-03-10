
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
