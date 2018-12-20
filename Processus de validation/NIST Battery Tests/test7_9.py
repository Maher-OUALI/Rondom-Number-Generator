#### teste et contient les conditions necessaires

import numpy as np
from scipy import special
from math import *
from itertools import *

#### 7eme test
#non overlapping Template Matching Test(igamc / chi square)
def nonoverlap_templ_match_test(ch,trgt='00000001',nbr_bloc=8): #trgt est un motif non periodic dont la frequence d'occurence est objet de ce test  
    """non-overlapping template match test: le but de ce test est de detecter la repitition d'un motif bien fixe """
    n=len(ch)
    N=nbr_bloc
    b=len(trgt)
    M=n//N
    var=M*((1/float(2**b))-((2*b-1)/float(2**(2*b))))
    moy=(M-b+1)/float(2**b)
    print(moy,var)#

    if(N<=100)and(b in range(2,11))and(M>0.01*n):
        L=[]
        for i in range(N):
            ch1=''
            for j in range(M):
                ch1=ch1+ch[i*M+j]
            L.append(ch1)
        L1=[]
        for i in range(N):
            s=0
            j=0
            while(j<=M):
                if(L[i][j:j+b]==trgt):
                    s+=1
                    j+=3
                else:
                    j+=1
            L1.append(s)
        Xobs=0
        for i in range(N):
            Xobs+=((L1[i]-moy)**2)/float(var)
        print('Xobs : ',Xobs)
        a=1-(special.gammainc(N/float(2),Xobs/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du septieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la sequence est de longueur insuffisante ')


#### 8eme test
#overlapping Template matching test(igamc / chi square distribution )
def overlap_templ_match_test(ch,trgt='11111111',long_bloc=1032,nbr_bloc=968):#trgt est un motif constitue que de 1 de longueur b
    """ overlapping template match test: le but de ce test est de detecter la repition d'un motif bien determine """
    n=len(ch)
    val_theor=[0.364091,0.185659,0.139381,0.100571,0.070432,0.139865]
    K=5#need to be 5
    M=long_bloc
    N=nbr_bloc
    b=len(trgt)
    lamda=(M-b+1)/float(2**b)
    nu=lamda*0.5
    if(N*0.070432>5)and(n>=10**4): 
        L=[]
        for i in range(N):
            ch1=''
            for j in range(M):
                ch1=ch1+ch[i*M+j]
            L.append(ch1)
        L1=[]
        for i in range(N):
            s=0
            j=0
            while(j<=M):
                if(L[i][j:j+b]==trgt):
                    s+=1
                j+=1
            L1.append(s)
        T=[0]*6
        for i in range(N):
            if(L1[i]>=5):
                T[5]+=1
            else:
                T[L1[i]]+=1
        Xobs=0
        for i in range(K+1):
            Xobs+=((T[i]-N*val_theor[i])**2)/float(N*val_theor[i])
        print('Xobs : ',Xobs)
        a=1-(special.gammainc(K/float(2),Xobs/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du huitieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la sequence est de longueur insuffisante')

    

#### 9eme test
#maurer test (erfc / half normal distribution )
def maurer_test(ch):
    
    #sous-fonctions utiles
    def toutes_combinaisons(m):
        """retourne toutes les sequences possibles de longueur m formees par 1 et 0 """
        lis=[]
        for i in product(*[[0,1]]*m):
            ch1=''
            for j in range(len(i)):
                ch1+=str(i[j])
            lis.append(ch1)
        return(lis)
    def pos_motif(L,mot):#L contient des motifs deux a deux distincts et mot existe presque surement dans L donc unicite et existence 
        """retourne la position d'un motif dans la liste des motifs """
        pos=None
        i=0
        while(i<=len(L))and(pos==None):
            if(L[i]==mot):
                pos=i
            i+=1
        return(pos)
    """ Murer Universal test: le but de ce test est de determiner """
    n=len(ch)
    val_theor=[(6,5.2177052,2.954),(7,6.1962507,3.125),(8,7.1836656,3.238),(9,8.1764248,3.311),(10,9.1723243,3.356),(11,10.170032,3.384),(12,11.168765,3.401),(13,12.168070,3.410),(14,13.167693,3.416),(15,14.167488,3.419),(16,15.167379,3.421)] #(L,expectedValue,Variance(L))
    L=6
    Q=10*(2**L) #voir les valeurs minimales que n peut prendre 
    K=n//L-Q
    c=0.7-0.8/float(L)+(4+32/float(L))*(K**(-3/float(L))/float(15))
    sig=c*sqrt(val_theor[L-6][2]/float(K))

    if(n>=101*Q*L)and(6<=L)and(L<=16):
        Lis=[]
        for i in range(n//L):
            ch1=''
            for j in range(L):
                ch1=ch1+ch[i*L+j]
            Lis.append(ch1)

        T=toutes_combinaisons(L)
        T1=[]
        for i in range(len(T)):
            s=0
            for j in range(Q):
                if T[i]==Lis[j]:
                    s=j+1
            T1.append((T[i],s))
        
        somme=0
        for i in range(Q,K+Q):
            pos=pos_motif(T,Lis[i])
            somme+=log(i+1-T1[pos][1],2)
            T1[pos]=(T[pos],i+1)
        somme=(somme/float(K))
        Vobs=abs((somme-val_theor[L-6][1])/float(sqrt(2*sig)))
        print('Vobs : ',Vobs)
        a=special.erfc(Vobs)
        res=bool(a>0.01)
        return(a,' / le resultat du neuvieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la sequence de longueur insuffisante ')
    
        
                
