#### teste et contient les conditions necessaires

from scipy import special
from math import *


##### 1er test
#monobit test(erfc // sur la sequence entiere // distribution normale  )
def monobittest(ch):
    """monobit test: determiner si la frequence des 1 est egale a celle des 0 """
    n=len(ch)
    if(n>=100):
        s=0
        for i in range(n):
            s=s+(2*int(ch[i])-1)
        k=(abs(s)/sqrt(n))
        print('abs(S)/sqrt(n) : ',k)
        a=(special.erfc(k/sqrt(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du premier test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la longueur de la sequence est insuffisante ')


#### 2eme test
#frequency within a block test(igamc  // sur M blocks of bits // distribution chi square )
def frequencybleocktest(ch,M=20):
    """frequency within a block test:determiner si la frequence des 1 est egale a celle des 0 dans des sous-sequences disjointes"""
    L=[]
    n=len(ch)
    p=n//M
    if(n>=M*p)and(n>=100)and(M>0.01*n)and(M>=20)and(p<100):
        for i in range(p):
            s=0
            for j in range(M):
                s+=int(ch[i*M+j])
            L.append(s/float(M))
        s=0
        for i in range(p):
            s+=(L[i]-0.5)**2
        K=4*M*s
        print('chi squared , ',K)
        a=1-(special.gammainc(p/float(2),K/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du deuxieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la longueur de la sequnece est insufisante')


#### 3eme test
#runs test(erfc // sur la sequence entiere // distribution chi square )
def runstest(ch):
    """runs test:caracterise l'oscillation entre des 0 et 1 en comptant le nombre des series non interrompues des 1 """
    n=len(ch)
    if(n>=100):
        s=0
        for i in range(n):
            s=s+int(ch[i])
        s=s/float(n)
        if(abs(s-0.5)<(2/sqrt(n))):
           v=1
           for i in range(n-1):
               v+=1-int(ch[i]==ch[i+1])
           print('Vn(obs) : ',v)
           k=abs(v-2*n*s*(1-s))
           k1=2*sqrt(2*n)*s*(1-s)
           a=special.erfc(k/float(k1))
           res=bool(a>0.01)
        else:
            a=0.00
            res=False
        return(a,' / le resultat du troisieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! la longueur de la sequence est insuffisante')
