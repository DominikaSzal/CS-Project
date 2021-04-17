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
import itertools as it
def assemble_genome(dna_list):
    
    for A in dna_list:
        
        if len(A)==1:
            return A[0]
        dic={}
        for i in range(len(A)):
            for j in range(len(A)):
                if i!=j:
                    ol=0
                    for k in range(1,min(len(A[i]),len(A[j]))):
                        if A[j][:k]==A[i][-k:]:
                            ol=k
                    dic[(i,j)]=ol
        if max(dic.values())==0:
            return "".join(A)
        else:
            ret="".join(A)
            l=len(ret)
            stack=[]
            for i,wd in enumerate(A):
                tmp=set(range(len(A)))
                tmp.remove(i)
                stack.append((wd,i,tmp))
            while stack:
                ans,cur,remain=stack.pop()
                if len(ans)<l:
                    if not remain:
                        ret=ans
                        l=len(ret)
                    else:
                        tmp=[[dic[cur,idx],idx] for idx in remain] # [#overlap,idx]
                        tmp.sort()
                        for ol,idx in tmp:
                            nans=ans+A[idx][ol:]
                            nremain=set(remain)
                            nremain.remove(idx)
                            stack.append((nans,idx,nremain))
            return ret
