def find_splice(dna_motif, dna):

    start = 0
    pos_list = []

    for char in dna_motif:
        position = dna.find(char,start)
  

        pos_list.append(position)
  

        start = position+1
        if position ==-1:
            return[]