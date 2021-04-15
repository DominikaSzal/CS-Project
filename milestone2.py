def find_splice(dna_motif, dna):

    start = 0
    pos_list = []

    for char in range( dna_motif):
       {
        position = dna.find(char,start)
  
        pos_list.append(position)
      
        position = 0
        while pos_list >15:
            position =dna.find(char,15)
            if char <=15 :
                return pos_list
        if dna_motif = False:
            return[]
         }   
    return pos_list
