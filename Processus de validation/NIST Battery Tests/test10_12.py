#### teste et comporte les conditions et valeurs limites 

import numpy as np
from scipy import special
from math import *
from itertools import *
from copy import *

##### 10eme test
#Linear complexity test(igamc / chi square distribution)

def complex_linea_test(ch,M=1000):
    #sous fonction utile
    def berlekamp_massey(ch1): #voir wikipedia
        a=len(ch1)
        b=[1]+[0]*(a-1)
        c=[1]+[0]*(a-1)
        l=0
        m=-1
        i=0
        L=[int(j) for j in ch1]
        while(i<a):
           v = L[(i - l):i]
           v = v[::-1]
           cc = c[1:l + 1]
           d = (L[i] + np.dot(v, cc)) % 2
           if d == 1:
              temp = copy(c)
              p = np.zeros(a)
              for j in range(0, l):
                if b[j] == 1:
                   p[j + i - m] = 1
              c = (c + p) % 2
              if l <= 0.5 * i:
                  l = i + 1 - l
                  m = i
                  b = temp
                   
           i+=1    
        return(l)
    n=len(ch)    
    val_theor=[0.010417,0.03125,0.125,0.5,0.25,0.0625,0.020833] # valeurs theorique 
    K=6 # fixed DOF
    N=n//M
    moy_theor=0.5*M+(9+(-1)**(M+1))/float(36)-(M/float(3)+2/float(9))/float(2**M)

    if(n>=10**6)and(M in range(500,5001))and(N<=200):
        L=[]
        for i in range(N):
            ch1=''
            for j in range(M):
                ch1=ch1+ch[i*M+j]
            L.append(ch1)
        

        lis=[]
        for i in range(N):
            lis.append(berlekamp_massey(L[i]))    
        

        T=[(-1)**M*(lis[i]-moy_theor)+2/float(9) for i in range(N)]
        fre=[0]*7 # de longueur K+1 => fixe 
        for i in range(N):
            if(T[i]<=-2.5):
                fre[o]+=1
            elif(T[i]<=-1.5):
                fre[1]+=1
            elif(T[i]<=-0.5):
                fre[2]+=1
            elif(T[i]<=0.5):
                fre[3]+=1
            elif(T[i]<=1.5):
                fre[4]+=1
            elif(T[i]<=2.5):
                fre[5]+=1
            else:
                fre[6]+=1
        Xobs=0
        for i in range(K+1):
            Xobs+=(fre[i]-N*val_theor[i])**2/N*val_theor[i]
        print('Xobs : ', Xobs)
        a=1-(special.gammainc(K/float(2),Xobs/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du dixieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! les approximations ne sont pas respectees ')

#### 11eme test
#Serial Test(igamc / chi square)
def serialTest(ch,m):#pour m=1 ce test est le meme que le 1er test #condition sur m,n m< int(log2(n))-2.
    #### des sous-fonctions utiles 
    def toutes_combinaisons(s):
        """retourne toutes les sequences possibles de longueur m formees par 1 et 0 """
        
        lis=[]
        for i in product(*[[0,1]]*s):
            ch1=''
            for j in range(len(i)):
                ch1+=str(i[j])
            lis.append(ch1)
        return(lis)
    def occurence(mot,ch):
        """ retourne le nombre d'occurence d'un motif dans une chaine de caractere en sutant un indice a chaque examination"""
        nbr=0
        b=len(mot)
        n=len(ch)
        i=0
        while(i+b<=n):
            if(mot==ch[i:i+b]):
                nbr+=1
            i+=1
        return(nbr)
    n=len(ch)
    if(m<int(log(n,2))-2):
        L=[]
        for i in range(m):
            L.append(ch+ch[0:m-i-1])
        print(L)
        lis=[]
        for i in range(m):
            lis.append(toutes_combinaisons(m-i))
        print(lis)
        T=[]
        for i in range(len(L)):
            temp=[]
            for j in range(len(lis[i])):
                temp.append(occurence(lis[i][j],L[i]))
                print(temp)
            T.append(temp)
        print('liste des occurences: ',T)#
        for i in range(len(T)):
            s=0
            for j in range(len(T[i])):
                s+=T[i][j]**2
            T[i]=(2**(m-i))/float(n)*s-n
        print('liste des frequences: ',T)#
        Xobs1=T[0]-T[1]
        Xobs2=T[0]-2*T[1]+T[2]
        print('Xobs1 : ',Xobs1,'   /   Xobs2 : ',Xobs2)
        a1=1-(special.gammainc(2**(m-2),Xobs1))
        a2=1-(special.gammainc(2**(m-3),Xobs2))
        res=bool(a1>0.01)and bool(a2>0.01)
        return(a1,' / ',a2,' / le resultat du onzieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la sequence est de longueur insuffisante ')


### 12 eme test
#Approximate Entropy Test(igamc / chi square distribution)
def Approx_entrop_test(ch,m=3):
    #### des sous-fonctions utiles 
    def toutes_combinaisons(s):
        """retourne toutes les sequences possibles de longueur m formees par 1 et 0 """
        lis=[]
        for i in product(*[[0,1]]*s):
            ch1=''
            for j in range(len(i)):
                ch1+=str(i[j])
            lis.append(ch1)
        return(lis)
    def occurence(mot,ch):
        """ retourne le nombre d'occurence d'un motif dans une chaine de caractere en sutant un indice a chaque examination"""
        nbr=0
        b=len(mot)
        n=len(ch)
        i=0
        while(i+b<=n):
            if(mot==ch[i:i+b]):
                nbr+=1
            i+=1
        return(nbr)
    n=len(ch)
    if(n>=1000)and(m<int(log(n,2))-5):
        L=[]
        for i in range(2):
            L.append(ch+ch[0:m+i-1])
        lis=[]
        for i in range(2):
            lis.append(toutes_combinaisons(m+i))
        T=[]
        for i in range(2):
            temp=[]
            for j in range(len(lis[i])):
                temp.append(occurence(lis[i][j],L[i]))
            T.append(temp)
        for i in range(2):
            s=0
            for j in range(len(T[i])):
                if(T[i][j] != 0):
                    T[i][j]=T[i][j]/float(n)
                    s+=(T[i][j])*log(T[i][j])
            T[i]=s
        Xobs1=2*n*(log(2)-T[0]+T[1])
        print('Xobs1 : ',Xobs1)
        a1=1-(special.gammainc(2**(m-1),Xobs1/float(2)))
        res=bool(a1>0.01)
        return(a1,' / le resultat du deuxieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la sequence est de longueur insuffisante ')
    
    
    
    
