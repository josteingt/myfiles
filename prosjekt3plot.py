# -*- coding: latin-1 -*-
from numpy import *

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as mpl


filnavn = "positions.txt"
telle = open(filnavn, "r")

#telle linjer
i = 0
for line in telle:
    i = i + 1
data = line.split()
k = len(data)
telle.close()

plotdata = zeros([k,i])

lesfil = open(filnavn, "r")
j = 0
for line in lesfil:
    data = line.split()
    for s in range(0,k):
        plotdata[s,j] = data[s]
    j = j+1

lesfil.close()


mpl.close("all")
mpl.figure(1)
for s in range(0,k,2):
    mpl.plot(plotdata[s],plotdata[s+1])

mpl.title("solsystem")
mpl.xlabel("AU")
mpl.ylabel("AU")
mpl.title("solsystemet")


filnavn = "energy.txt"
telle = open(filnavn, "r")

#telle linjer
k = 0
for line in telle:
    k = k + 1

telle.close()

energy = zeros(k)

lesfil = open(filnavn, "r")

t = linspace(0,k,k)

i = 0
for line in lesfil:
    energy[i] = line
    i = i + 1
lesfil.close()

mpl.figure(50)
mpl.plot(t,energy)
mpl.title("total energy")
mpl.ylabel("energy")
mpl.xlabel("number of timesteps")

mpl.show()
