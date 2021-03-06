##################################################### Import system libraries ######################################################
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds)
from scipy.optimize import curve_fit


#Auslesen der Daten und Umrechnen in die Magnetfeldstärken
f, pk1, pk2, I_Horiz1,I_Horiz2  = np.genfromtxt('messdaten/aufg_c.txt', unpack=True)

#Umrechnen in die passenden Stromstärken

I_sweep_1, I_sweep_2= pk1*0.1,  pk2*0.1
I_Horiz1, I_Horiz2=I_Horiz1*0.3, I_Horiz2*0.3

#Daten der Spulen
#Sweep Spule

R_sweep, N_sweep= 0.1639, 11
R_hori, N_hori= 0.1579, 154


#Umrechnung der Stromstärken in das Magnetfeld in mT

def Bcalc(I,N, R):
    return 8*const.mu_0*I*N/(np.sqrt(125)*R)
B_peak1=(Bcalc(I_sweep_1,N_sweep, R_sweep)+Bcalc(I_Horiz1, N_hori, R_hori))
B_peak2=(Bcalc(I_sweep_2,N_sweep, R_sweep)+Bcalc(I_Horiz2, N_hori, R_hori))


#Die Funktionsanpassungen
def fitfunction(x, m, b):
    return m*x + b

params_1, covariance_1 = curve_fit( fitfunction, f, B_peak1 , maxfev=900000000)
params_2, covariance_2 = curve_fit( fitfunction, f, B_peak2 , maxfev=900000000)
errors_1 = np.sqrt(np.diag(covariance_1))
errors_2 = np.sqrt(np.diag(covariance_2))

#Print der Parameter
print("Parameter der Fits")
print("m_1=0", params_1[0], "±", errors_1[0])
print("b_1=0", params_1[1], "±", errors_1[1])
print("m_2=0", params_2[0], "±", errors_2[0])
print("b_1=0", params_2[1], "±", errors_2[1])


#Plotten
x_plot = np.linspace(0, 1000 ,10000)
plt.plot(x_plot, fitfunction(x_plot, *params_1), "r-", label="Fit_1")
plt.plot(x_plot, fitfunction(x_plot, *params_2), "r-", label="Fit_2")
plt.plot( f, B_peak1, "kx",label="Messwerte des ersten Peaks")
plt.plot( f, B_peak2, "kx",label="Messwerte des zweiten Peaks")
plt.xlabel('f in kHz')
plt.ylabel('B in T')
plt.legend(loc="best")
plt.savefig("ressources\fit_reso.pdf")


#Berechung der Lande-Faktoren
#Konstanten vom Bruder Scipy mit Unsicherheit, die aber vernachlässigbar ist
(value_mub,unit_mub,uncertainty)=const.physical_constants["Bohr magneton"]

h=const.h
mu_b=value_mub
g1, g2=h/(ufloat(params_1[0],errors_1[0])*mu_b), h/(ufloat(params_2[0],errors_2[0])*mu_b)
print("g1,g2")
print(g1,g2)
I_vert=ufloat(0.15,0.01)

#Berechnung des vertikalen Magnetfeldes
B_verti=I_vert*8*const.mu_0*20/(np.sqrt(125)*0.117)
B_hori1=ufloat(params_1[1],errors_1[1])
B_hori2=ufloat(params_1[1],errors_1[1])
B_hori=(B_hori1+B_hori2)/2
#Berechnung des gesamten Erdmagnetfeldes
B_ges=(B_hori**2+B_verti**2)
print("B_gesamt",B_ges)
print("B_hori",B_hori)
print("B_verti",B_verti)
