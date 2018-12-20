#### teste et comporte les conditions et valeurs limites 

import numpy as np
from scipy import special
from math import *
from itertools import *
from scipy.stats import norm #on va utiliser la norm.cdf pour calculer la " cumultative 


#### 13eme test
#Cumulative Sums (Cusum) Test(erfc / normal distribution )
def somme_cumul(ch,mode=0): #mode=0 => forward , mode=1 => backward    #teste
    n=len(ch)
    if(n>=100):
        L=[]# liste des sommes partielles
        s=0
        for i in range(n)[::-1*int(mode==1)+int(mode==0)]:
            s+=2*int(ch[i])-1
            L.append(s)
     
        L1=[abs(L[i]) for i in range(len(L))]
        z=max(L1)
        a=int((-n/float(z)+1)/float(4))
        b=int((-n/float(z)-3)/float(4))
        c=int((n/float(z)-1)/float(4))
     
        Xobs=0
        for i in range(a,c+1):
            Xobs+=norm.cdf((4*i+1)*z/float(sqrt(n)))-norm.cdf((4*i-1)*z/float(sqrt(n)))
        print(Xobs)
        for i in range(b,c+1):
            Xobs-=norm.cdf((4*i+3)*z/float(sqrt(n)))-norm.cdf((4*i+1)*z/float(sqrt(n)))
            
        a=1-Xobs
        res=bool(a>=0.01)
        return(a,' / le resultat du treizieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la sequence est de longueur insuffisante')


####14eme test
#Random Excursions Test(igamc / chi square distribution )
def rand_excursion_test(ch):
    ##sous-fonctions utiles
    def occurence_liste(L,x):
        res=[]
        for i in range(len(L)):
            s=0
            for j in range(len(L[i])):
                if(L[i][j]==x):
                   s+=1
            res.append(s)
        return(res)
    n=len(ch)
    val_theor=[[0.5,0.25,0.125,0.0625,0.0312,0.0312],[0.75,0.0625,0.0469,0.0352,0.0264,0.0791],[0.8333,0.0278,0.0231,0.0193,0.0161,0.0804],[0.8750,0.0156,0.0137,0.012,0.0105,0.0733],[0.9,0.01,0.009,0.0081,0.0073,0.0656],[0.9167,0.0069,0.0064,0.0058,0.0053,0.0588],[0.9286,0.0051,0.0047,0.0044,0.0041,0.531]]
    #ces valeurs sont calculees a partir d'une loi geometrique pour sept etat sachant que P(x)=P(-x)
    if(n>=10**6):
        L=[]
        for i in range(n):
            L.append(2*int(ch[i])-1)# list formee par des 1 et -1 a la place des 0
        L1=[] # liste des sommes partielles
        s=0
        for i in range(n):
            s+=L[i]
            L1.append(s)
        L1.insert(0,0)
        L1.append(0)  ##attention len(L1)=n+2
        L=[] # liste des cycles 
        temp=[]
        for i in range(n+2):
            temp.append(L1[i])
            if(L1[i]==0)and(i!=0): 
                L.append(temp)
                temp=[0]
        J=len(L)
        T=[]# liste des occurences des cycles 
        for i in range(4):
            temp=[]
            for j in range(len(L)):
                s=0
                for k in range(len(L[j])):
                    if(L[j][k]==i-4):
                        s+=1
                temp.append(s)
            T.append(temp)
        for i in range(4):
            temp=[]
            for j in range(len(L)):
                s=0
                for k in range(len(L[j])):
                    if(L[j][k]==i+1):
                        s+=1
                temp.append(s)
            T.append(temp)
        freq=[] #liste des frequence nu i 
        for i in range(4):
            temp=[0]*6
            for j in range(len(T[i])):
                if(T[i][j]>=5):
                    temp[5]+=1
                else:
                    temp[T[i][j]]+=1            
            freq.append(temp)
        for i in range(4):
            temp=[0]*6
            for j in range(len(T[i+4])):
                if(T[i+4][j]>=5):
                    temp[5]+=1
                else:
                    temp[T[i+4][j]]+=1            
            freq.append(temp)
        print('liste des frequences ',freq)#
        Xobs=[]
        for i in range(4):
            s=0
            for j in range(6):
                s+=(freq[i][j]-J*val_theor[3-i][j])**2/float(J*val_theor[3-i][j])
            Xobs.append(s)
        for i in range(4):
            s=0
            for j in range(6):
                s+=(freq[i+4][j]-J*val_theor[i][j])**2/float(J*val_theor[i][j])    ######need to be checked
            Xobs.append(s)
        print('liste des valeurs observes ',Xobs)#
        P=[special.gammainc(5/2.0,Xobs[i]/2.0) for i in range(len(Xobs))]
        res=1
        for i in range(8):
            res*=bool(P[i]>=0.01)
        
        return(P,' / le resultat du quinzieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la longueur de la sequence est insuffisante')
            
    
#### 15eme test
#Random Excursions Variant Test(erfc / normal distribution /done on the hole sequence )    #teste
def rand_excursion_variant_test(ch):
    #des sous-fonctions utiles 
    def occurence_(L,x):
        """retourner le nombre d'occurence d'un certin element dans un tableau"""
        s=0
        for i in range(len(L)):
            if L[i]==x:
                s+=1
        return(s)
    n=len(ch)
    if(n>=10**6):
        L=[]
        for i in range(n):
            L.append(2*int(ch[i])-1)# list formee par des 1 et -1 a la place des 0
        L1=[] # liste des sommes partielles
        s=0
        for i in range(n):
            s+=L[i]
            L1.append(s)
        L1.insert(0,0)
        L1.append(0)  ##attention len(L1)=n+2
        J=0
        for i in range(1,n+2):
            if(L1[i]==0):
                J+=1
        T=[(i-9,occurence_(L1,i-9)) for i in range(9)]
        T=T+[(i+1,occurence_(L1,i+1)) for i in range(9)]
        T1=[abs((T[i][1]-J)/float(sqrt(2*J*(4*abs(T[i][0])-2)))) for i in range(18)]
        T1=special.erfc(T1)
        T1=list(T1)
        res=1
        for i in range(18):
            res*=bool(T1[i]>=0.01)
        
        return(T1,' / le resultat du quinzieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la longueur de la sequence est insuffisante ')
        
        
    
    
    
