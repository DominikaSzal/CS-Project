def find_splice(dna_motif, dna):

    start = 0
    pos_list = []

    for char in dna_motif:
        position = dna.find(char,start)
  

        pos_list.append(position)
  

        start = position+1
        if position ==-1:
            return[]
    return pos_list
def assemble_genome(A):

    N = len(A)

 

    # Populate overlaps

    overlaps = [[0] * N for _ in range(N)]

    for i, x in enumerate(A):

        for j, y in enumerate(A):

            if i != j:

                for ans in range(min(len(x), len(y)), -1, -1):

                    if x.endswith(y[:ans]):

                        overlaps[i][j] = ans

                        break

 

    # dp[mask][i] = most overlap with mask, ending with ith element

    dp = [[0] * N for _ in range(1<<N)]

    parent = [[None] * N for _ in range(1<<N)]

    for mask in range(1, 1 << N):

        for bit in range(N):

            if (mask >> bit) & 1:

                # Let's try to find dp[mask][bit].  Previously, we had

                # a collection of items represented by pmask.

                pmask = mask ^ (1 << bit)

                if pmask == 0: continue

                for i in range(N):

                    if (pmask >> i) & 1:

                        # For each bit i in pmask, calculate the value

                        # if we ended with word i, then added word 'bit'.

                        value = dp[pmask][i] + overlaps[i][bit]

                        if value > dp[mask][bit]:

                            dp[mask][bit] = value

                            parent[mask][bit] = i

 

    # Answer will have length sum(len(A[i]) for i) - max(dp[-1])

    # Reconstruct answer:

 

    # Follow parents down backwards path that retains maximum overlap

    perm = []

    mask = (1<<N) - 1

    i = max(range(N), key = dp[-1].__getitem__)

    while i is not None:

        perm.append(i)

        mask, i = mask ^ (1<<i), parent[mask][i]

 

    # Reverse path to get forwards direction; add all remaining words

    perm = perm[::-1]

    seen = [False] * N

    for x in perm:

        seen[x] = True

    perm.extend([i for i in range(N) if not seen[i]])

 

    # Reconstruct answer given perm = word indices in left to right order

    ans = [A[perm[0]]]

    for i in range(1, len(perm)):

        overlap = overlaps[perm[i-1]][perm[i]]

        ans.append(A[perm[i]][overlap:])

 

    return "".join(ans)

def shared_motif(dna_list):
    total = len(dna_list)
    #using first string as a reference
    ref = dna_list[0]
    ref_len = len(ref)

    substring = '' 
    
    if total > 1 and ref_len > 0:
       for i in range(ref_len):
          for d in range((ref_len)-i+1):
             if d > len(substring) and all(ref[i:i+d] in x for x in dna_list):
                substring = ref[i:i+d]



    return substring 


def reverse_complement(dna):
    dna = dna[::-1]
    reverse_dna = ''
    
    for symbol in dna:
        if symbol == 'A':
            reverse_dna +='T'
        elif symbol == 'C':
            reverse_dna += 'G'
        elif symbol == 'G':
            reverse_dna +='C'
        elif symbol == 'T':
            reverse_dna +='A'
    return reverse_dna
def rev_palindrome(dna):
    ans = []

  
    for i in range(len(dna)-4):
        for j in range(i+3,min(len(dna),i+12)):
            
            string = dna[i:j+1]
            if reverse_complement(dna[i:j+1]) == string:
                ans.append((i,j-i+1))

    return ans

import math   
def random_genome(dna, gc_content):
   dna = dna.upper()
   GC = len(dna.replace('A', '').replace('T', ''))
   AT = len(dna.replace('C', '').replace('G', ''))
   #creating empty list to put logs in
   answer = []
   for i in range(0, len(gc_content)):
       x = GC * math.log10(float(gc_content[i]) / 2) + AT * math.log10((1 - float(gc_content[i])) / 2)
       answer.append(x)
   return answer


import numpy as np
from collections import Counter
import math

def perfect_match(rna):
    counts = Counter(rna)
    if counts['A'] == counts['U'] and counts['G'] == counts['C']:
        return math.factorial(counts['A'])*math.factorial(counts['A'])

def get_edgs(dict):
        #Initializing Variables
        Ros_IDs = dict.keys()
        RI_list=[]
        answer=[]

        #Creating a list of the ROSALIND identifiers
        for i in Ros_IDs:
            RI_list.append(i)

        #Compare the 3 figure prefix/suffix of items in RI_list to see if they share a common edge
        for i in range(0,len(RI_list)):
            for j in range(i+1,len(RI_list)):
                if(dict[RI_list[i]][:3]==dict[RI_list[j]][-3:] or dict[RI_list[i]][-3:]==dict[RI_list[j]][:3]):

                    #Add adjacentcy to the return list
                    answer.append((RI_list[i],RI_list[j]))   

        return answer
    
