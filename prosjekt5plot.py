# -*- coding: latin-1 -*-
from numpy import *

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as mpl

tid = "10"

filnavn = "eksplisitt.txt"
telle = open(filnavn, "r")

#telle linjer
i = 0
for line in telle:
    i = i + 1
telle.close()

posisjonex = zeros(i)
funksjonex = zeros(i)

lesfil = open(filnavn, "r")
j = 0
for line in lesfil:
    [pos,funk] = line.split()
    posisjonex[j] = pos
    funksjonex[j]= funk
    j = j+1

posisjonimp = zeros(i)
funksjonimp = zeros(i)

filnavn = "implisitt.txt"
lesfil = open(filnavn, "r")
j = 0
for line in lesfil:
    [pos,funk] = line.split()
    posisjonimp[j] = pos
    funksjonimp[j]= funk
    j = j+1

posisjoncn = zeros(i)
funksjoncn = zeros(i)

filnavn = "cn.txt"
lesfil = open(filnavn, "r")
j = 0
for line in lesfil:
    [pos,funk] = line.split()
    posisjoncn[j] = pos
    funksjoncn[j]= funk
    j = j+1

posisjonana = zeros(i)
funksjonana = zeros(i)

filnavn = "analytical.txt"
lesfil = open(filnavn, "r")
j = 0
for line in lesfil:
    [pos,funk] = line.split()
    posisjonana[j] = pos
    funksjonana[j]= funk
    j = j+1


lesfil.close()
mpl.figure(1)   
mpl.plot(posisjonex,funksjonex, "r")
mpl.plot(posisjonimp, funksjonimp, "k")
mpl.plot(posisjoncn, funksjoncn, "c")
mpl.plot(posisjonana, funksjonana,"y")
mpl.legend(["eksplisitt", "implisitt", "crank-nicolson", "analytisk"])
mpl.xlabel("posisjon")
mpl.ylabel("funksjonsverdi")
mpl.title("u(x,t=10)")
#mpl.savefig('omega' + omega + "without.png") 

residualex = zeros(i)
residualimp = zeros(i)
residualcn = zeros(i)

for i in range(i):
    residualex[i] = funksjonex[i] - funksjonana[i]
    residualimp[i] = funksjonimp[i] - funksjonana[i]
    residualcn[i] = funksjoncn[i] - funksjonana[i]

mpl.figure(2)
mpl.plot(posisjonex, residualex)
mpl.plot(posisjonimp, residualimp)
mpl.plot(posisjoncn, residualcn)
mpl.legend(["eksplisitt","implisitt", "crank-nicolson"])
mpl.xlabel("posisjon")
mpl.ylabel("funksjonsverdi")
mpl.title("residualplot")
mpl.show()
