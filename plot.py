# -*- coding: latin-1 -*-
from numpy import *

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as mpl


omega = "0_01"

filnavn = "psikvadrat"+omega+"without.txt"
telle = open(filnavn, "r")

#telle linjer
i = 0
for line in telle:
    i = i + 1
telle.close()

psikvadrat = zeros(i)
posisjon = zeros(i)

lesfil = open(filnavn, "r")
j = 0
for line in lesfil:
    [p,psi] = line.split()
    posisjon[j]= p
    psikvadrat[j] = psi
    j = j+1

lesfil.close()
   
mpl.plot(posisjon,psikvadrat)
mpl.xlabel("r/a")
mpl.ylabel("psikvadrat")
mpl.title("omega=" + omega)
mpl.savefig("psikvadrat" + omega + "without.png")
mpl.show()
