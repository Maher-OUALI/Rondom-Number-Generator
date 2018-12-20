from PIL import Image
import numpy as np
import random as rd
from math import sqrt,cos,sin,pi
import matplotlib.pyplot as plt

ch=''
for i in range(10000):
        ch+=str(rd.randint(0,1))
        
        
### test visuel bitmap image d'une sequence aleatoire 
def bitmap(ch):
        n=len(ch)
        h=int(sqrt(n))
        data = np.zeros((h, h, 3), dtype=np.uint8)
        for i in range(h):
                for j in range(h):
                        data[i,j]=[255*(1-int(ch[h*i+j]))]*3


        img = Image.fromarray(data, 'RGB')
        img.save('my.png')
        img.show()

### test approximatif de la valeur de pi d'une sequence aleatoire selon la methode de monte carlo
def circle():
        X,Y=[],[]
        for i in range(1000):
                X.append(cos(pi*0.5*(i/float(1000))))
                Y.append(sin(pi*0.5*(i/float(1000))))
        plt.plot(X,Y,'blue')
def montecarlo(ch):
        p=len(ch)//8
        T=[]
        for i in range(p):
                if(i%2==0):
                        a1=0
                        for j in range(8):
                                a1+=int(ch[8*i+j])*(1/float(2**(j+1)))
                else:
                        a2=0
                        for j in range(8):
                                a2+=int(ch[8*i+j])*(1/float(2**(j+1)))
                        T.append((a1,a2))
        somme_out,somme_in=0,0
        X_in,Y_in,X_out,Y_out=[],[],[],[]
        for i in range(len(T)):
                if(T[i][0]**2+T[i][1]**2>1):
                        somme_out+=1
                        X_out.append(T[i][0])
                        Y_out.append(T[i][1])
                else:
                        somme_in+=1
                        X_in.append(T[i][0])
                        Y_in.append(T[i][1])
        plt.plot(X_in,Y_in,'go')
        plt.plot(X_out,Y_out,'ro')
        circle()
        plt.show()
        n=len(T)
        #somme_in/n est equivalente a surface_cercle/surface_carre = (1/4 * pi *R**2)/R**2 = 1/4 * pi => pi = 4 *somme_in /n
        approx_pi=4*somme_in/float(n)
        print('la valeur approximative de pi selon la methode monte carlo est ',approx_pi)

### programme principal               
print(ch)
bitmap(ch)
montecarlo(ch)

