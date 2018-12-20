#### teste et contient les conditions necessaires

from scipy import special
from math import *
import numpy as np
from scipy import fftpack as f

#####4eme test 
#longest run of ones within a block test(8) (igamc  // sur M blocks of bits // distribution chi square )
def longest_runs_of_ones_8(ch):
    """longest runs of ones within 8_bits blocks: caracterise l'oscillation entre des 0 et des 1 dans des sous-sequences de longueur 8"""
    val_theor=[0.2148, 0.3672, 0.2305, 0.1875]#valeurs theoriques calculees dans le cas d'une sequence aleatoires
    n=len(ch)
    k,M=3,8      #valeurs theoriques
    if(n>=128): #valeur minimale de n 
        L=[]
        p=n//M
        for i in range(p):
            ch1=''
            for j in range(M):
                ch1=ch1+ch[i*M+j]
            L.append(ch1)    #division en M_bits blocs
        
        L1=[]
        for i in range(p):
            m=0
            s=0
            for j in range(M):
                if (L[i][j]=='1'):
                    s=s+1
                else:
                     s=0
                if s>m:    
                    m=s
            L1.append(m) #calcul de la sequence de 1 successives la plus longue dans chaque bloc
             
        L_freq=[0,0,0,0]
        for i in range(p):
            if (L1[i]<=1):
                L_freq[0]+=1
            elif(L1[i]>=4):
                L_freq[3]+=1
            else:
                L_freq[L1[i]-1]+=1
        Xobs=0
        for i in range(k+1):
            Xobs+=(L_freq[i]-p*val_theor[i])**2/float(p*val_theor[i])
        print('Xobs : ',Xobs)
        a=1-(special.gammainc(k/float(2),Xobs/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du quatrieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !!! la longueur de la sequence est insuffisant')

#longest run of ones within a block test(128) (igamc  // sur M blocks of bits // distribution chi square )
def longest_runs_of_ones_128(ch,n):
    """longest runs of ones within 128_bits blocks:caracterise l'oscillation entre des 0 et des 1 dans des sous-sequences de longueur 128"""
    val_theor=[ 0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124 ]#valeurs theoriques calculees dans le cas d'une sequence aleatoires
    M,k=128,5       #valeurs theoriques

    if(n>=6272):
        L=[]
        p=n//M
        for i in range(p):
            ch1=''
            for j in range(M):
                ch1=ch1+ch[i*M+j]
            L.append(ch1)    #division en M_bits blocs
       
        L1=[]
        for i in range(p):
            m=0
            s=0
            for j in range(M):
                if (L[i][j]=='1'):
                    s=s+1
                else:
                     s=0
                if s>m:    
                    m=s
            L1.append(m) #calcul de la sequence de 1 successives la plus longue dans chaque bloc
           
        L_freq=[0,0,0,0,0,0]
        for i in range(p):
            if (L1[i]<=4):
                L_freq[0]+=1
            elif(L1[i]>=9):
                L_freq[5]+=1
            else:
                L_freq[L1[i]-4]+=1
        Xobs=0
        for i in range(k+1):
            Xobs+=(L_freq[i]-p*val_theor[i])**2/float(p*val_theor[i])
        print('Xobs : ',Xobs)
        a=1-(special.gammainc(k/float(2),Xobs/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du quatrieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !!! la longueur de la sequence est insuffisant')
    
#longest run of ones within a block test(10000) (igamc  // sur M blocks of bits // distribution chi square )
def longest_runs_of_ones_10000(ch):
    """longest runs of ones within 10000_bits blocks: caracterise l'oscillation entre des 0 et des 1 dans des sous-sequences de longueur 10000"""
    val_theor=[0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]#valeurs theoriques calculees dans le cas d'une sequence aleatoires
    M,k=10000,6         #valeurs theoriques
    n=len(ch)

    if(n>=750,000):
        L=[]
        p=n//M
        for i in range(p):
            ch1=''
            for j in range(M):
                ch1=ch1+ch[i*M+j]
            L.append(ch1)    #division en M_bits blocs

        L1=[]
        for i in range(p):
            m=0
            s=0
            for j in range(M):
                if (L[i][j]=='1'):
                    s=s+1
                else:
                     s=0
                if s>m:    
                    m=s
            L1.append(m) #calcul de la sequence de 1 successives la plus longue dans chaque bloc
            
        L_freq=[0,0,0,0,0,0,0]
        for i in range(p):
            if (L1[i]<=10):
                L_freq[0]+=1
            elif(L1[i]>=16):
                L_freq[6]+=1
            else:
                L_freq[L1[i]-10]+=1
        Xobs=0
        for i in range(k+1):
            Xobs+=(L_freq[i]-p*val_theor[i])**2/float(p*val_theor[i])
        print('Xobs : ',Xobs)
        a=1-(special.gammainc(k/float(2),Xobs/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du quatrieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !!! la longueur de la sequence est insuffisant')


##### 5eme  test 
#Binary Matrix test(igamc /  chi square distribution)
def matrice_binaire_test(ch):
    """Binary Matrix Test : ce test detecte la dependance lineaire entre les sous-sequences de meme longueur fixe M"""
    val_theor=[0.2888,0.5776,0.1336] #
    n=len(ch)
    Q,M=32,32 #ligne
     #colonne   #valeurs theoriques due aux approximations considerees
    if(n>=38*M*Q):
        L=[]
        p=n//(M*Q)
        
        l=0
        for i in range(p):
            L1=[]
            for j in range(M):
                L2=[]
                for k in range(Q):
                    L2.append(int(ch[l]))
                    l+=1
                L1.append(L2)
            L.append(np.matrix(L1))  
        L2=[np.linalg.matrix_rank(L[i]) for i in range(len(L))]
        L1=[0,0,0]
        for i in range(p):
            if (L2[i]==M):
                L1[0]+=1
            elif(L2[i]==M-1):
                L1[1]+=1
            else:
                L1[2]+=1
        Xobs=((L1[0]-p*val_theor[0])**2/float(p*val_theor[0]))+((L1[1]-p*val_theor[1])**2/float(p*val_theor[1]))+((L1[2]-p*val_theor[2])**2/float(p*val_theor[2]))
        print('Xobs : ',Xobs)
        a=1-(special.gammainc(1,Xobs/float(2)))
        res=bool(a>0.01)
        return(a,' / le resultat du cinqueme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! longueur de sequence insuffisante ')



##### 6eme test
#Discrete Fourier transfrom test(erfc / normal distribution )
def FFTtest(ch):
    """Fast Fourier Transform test: determine si le nombre d'amplitude des frequence de la TF de la sequence qui depassse 95% d'une valeur reference est different de 5%"""
    n=len(ch)
    if(n>=1000):
        thres=sqrt(-log(0.05)*n)
        N0=0.95*0.5*n

        L=[2*int(ch[i])-1 for i in range(n)]
        A=f.fft(L)
        A=abs(A)[1:(n/2):] #need more work or else ghamma
        N1=0
        for i in range(len(A)):
            if (A[i]<thres):
                N1+=1
        d=(N1-N0)/float(sqrt(n*0.25*0.95*0.05))
        print('d : ',d)
        a=(special.erfc(abs(d)/float(sqrt(2))))
        res=bool(a>0.01)
        return(a,' / le resultat du sixieme test est ',int(res)*'random'+(1-int(res))*'non_random')
    else:
        return('ERROR !! longueur de sequence insuffisante pour avoir de vraies resultats ')

