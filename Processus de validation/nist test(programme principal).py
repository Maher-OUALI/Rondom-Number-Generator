#importation des modules et des fichiers des tests
import random
from scipy import special
from math import *
from test1_3 import *
from test4_6 import *
from test7_9 import *
from test10_12 import *
from test13_15 import *

#generateur de sequence aleatoire avec la commande random.randint(0,1)
ch=''
n=1000000  #longueur de la sequence de bits aleatoires
for i in range(n):
    x=random.randint(0,1)
    ch+=str(x)


print(ch)#affichage de la sequence 

#affichage des resultats des tests 
print(monobittest(ch))
print(frequencybleocktest(ch))
print(runstest(ch))
print(longest_runs_of_ones_128(ch))
print(matrice_binaire_test(ch))
print(FFTtest(ch))
print(nonoverlap_templ_match_test(ch,'00001'))
print(overlap_templ_match_test(ch,'11111'))
print(maurer_test(ch))
print()
print(serialTest(ch,3))
print(Approx_entrop_test(ch,3))
print(somme_cumul(ch,0))
print(rand_excursion_test(ch))
print(rand_excursion_variant_test(ch))
