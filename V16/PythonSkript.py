##################################################### Import system libraries ######################################################
import matplotlib as mpl
mpl.rcdefaults()
mpl.rcParams.update(mpl.rc_params_from_file('meine-matplotlibrc'))
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds,
)
################################################ Finish importing system libraries #################################################

################################################ Adding subfolder to system's path #################################################
import os
import sys
import inspect
# realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0], "python_custom_scripts")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
############################################# Finish adding subfolder to system's path #############################################

##################################################### Import custom libraries ######################################################
from curve_fit import ucurve_fit
from table import (
    make_table,
    make_full_table,
    make_composed_table,
    make_SI,
    write,
    search_replace_within_file,
)
from regression import (
    reg_linear,
    reg_quadratic,
    reg_cubic
)
from error_calculation import(
    mean,
    MeanError
)
from utility import(
    constant
)
################################################ Finish importing custom libraries #################################################
# epsilon = 8.854 * 10**(-12)
# e = 1.602 * 10**(-19)
# m_0 = 6.645 * 10**(-27)  # kg
# v = 16.3 * 10**6  # m/s
# I = 15.62  # eV
# R = 296.8  # J/kg k
# T = 293.15  # K
# z = 2
# Z = 7

# ln = np.log(2 * m_0 * v**2 / I)


# def Druck(p):
#     return 4 * np.pi * m_0 * v**2 * epsilon**2 * R * T * 100 / (e**4 * z**2 * Z * p * 10**5 * ln)


# p_plot = np.linspace(0, 1, 100)
# plt.plot(p_plot, Druck(p_plot), 'r-', label='Bla')
# plt.xlabel("p / [mbar]")
# plt.ylabel("U / [V]")
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/Druck.pdf')
# plt.clf()

#############################################

e = 1.602176 * 10**(-19)  # C
Z = 7
z = 2
m0J = 0.511 * 1.602 * 10**(-13)  # J
m0eV = 0.511 * 10**6  # eV
c = 299792458
e0 = 8.8541878 * 10**(-12)  # epsilon 0
v = 16.3 * 10**6
vstrich = 0.05425  # vstrich=v/c
I = 15.62  # eV
pi = 3.141592654
Ealpha = 5.486 * 1.602 * 10**(-13)  # J
NA = 6.02214178 * 10**23  # 1/mol
MLuft = 14 * 10**(-3)  # kg/mol Molare Masse (78% N2, 21% O2, 1% vernachlÃ¤ssigt)
RS = 296.8  # J/(kg K) spezifische Gaskonstante


def R(p):  # Druck in mbar --> R in cm
    return m0J * vstrich**2 * e0**2 * 4 * pi / (e**4 * z**2 * Z) * 1 / (np.log(2 * m0eV * vstrich**2 / I)) * Ealpha * 1 / ((NA * p) / (MLuft * RS * 293.15))


write('build/Reichweite_Luft.tex', make_SI(R(1013.25), r'\centi\meter', figures=1))
print(R(1013.25))
p = np.linspace(1, 1013.15, 1000)
plt.plot(p, R(p), 'b-', label='Reichweite')
plt.xlabel('Druck / [mbar]')
plt.ylabel('Reichweite / [cm]')
plt.tight_layout()
plt.legend(loc='best')
plt.xlim(0, 1014)
plt.ylim(0, 100)
# plt.grid()
plt.savefig('build/Druck.pdf')
plt.clf()


#############################################
Akt = 5608 / 300
write('build/Akt.tex', make_SI(Akt, r'\becquerel', figures=1))


Druck, U_max, U_min = np.genfromtxt('messdaten/ohne_folie.txt', unpack=True)
Druck_folie, U_folie_max, U_folie_min = np.genfromtxt('messdaten/mit_folie.txt', unpack=True)


U_mean = []
U_folie_mean = []
# for i in range(len(U_mean)):
#     U_mean[i] = np.mean(U_max[i], U_min[i])

for i in range(len(U_max)):
    a = []
    a.append(U_max[i])
    a.append(U_min[i])
    U_mean.append(a)

for i in range(len(U_folie_max)):
    a = []
    a.append(U_folie_max[i])
    a.append(U_folie_min[i])
    U_folie_mean.append(a)


U = (mean(U_mean, axis=1))
U_folie = (mean(U_folie_mean, axis=1))


t_plot = np.linspace(0, 280, 2)
params = ucurve_fit(reg_linear, Druck, noms(U))
params1 = ucurve_fit(reg_linear, Druck_folie, noms(U_folie))
write('build/a_ohne_Folie.tex', make_SI(params[0], r'\volt\per\milli\bar', figures=1))
write('build/a_mit_Folie.tex', make_SI(params1[0], r'\volt\per\milli\bar', figures=1))
write('build/b_ohne_Folie.tex', make_SI(params[1], r'\volt', figures=1))
write('build/b_mit_Folie.tex', make_SI(params1[1], r'\volt', figures=1))

plt.plot(t_plot, reg_linear(t_plot, *noms(params)), 'b-', label='Fit(ohne Folie)')
plt.errorbar(Druck, noms(U), fmt='bx', yerr=stds(U), label='Messdaten ohne Folie')        # mit Fehlerbalken

plt.plot(t_plot, reg_linear(t_plot, *noms(params1)), 'r-', label='Fit(mit Folie)')
plt.errorbar(Druck_folie, noms(U_folie), fmt='rx', yerr=stds(U_folie), label='Messdaten mit Folie')        # mit Fehlerbalken

plt.xlabel("p / [mbar]")
plt.ylabel("U / [V]")
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/test.pdf')
plt.clf()

print('HIER')
delta_E = 5.486 * 10**6 * 1.602 * 10**(-19) * (1 - (params1[1] / params[1]))  # in Joule
print(5.486 * 10**6 * (1 - (params1[1] / params[1])))
N = 5.9 * 10**28
Z = 79
I = 790  # eV
# v = 0.054  # c
v = unp.sqrt(5.486 * 10**6 * (1 + (params1[1] / params[1])) / (3727.38 * 10**6))
print(v)
z = 2
epsilon = 8.854 * 10**(-12)
e = 1.602 * 10**(-19)
m_0 = 511000  # eV/c^2


d = (m_0 * e * v**2 * 4 * np.pi * epsilon**2 * delta_E) / (e**4 * N * z**2 * Z * unp.log(2 * m_0 * v**2 / I))
write('build/Dicke_Folie.tex', make_SI(d * 10**(6), r'\micro\meter', figures=1))
print(d)

################################
# Wirkungsquerschnitt

Winkel, Counts, t = np.genfromtxt('messdaten/winkel.txt', unpack=True)

Counts_norm = Counts / t  # normiere die Counts
Counts_fehler = np.sqrt(Counts) / t
Counts = unp.uarray(Counts_norm, Counts_fehler)
print('Fehler')
print(Counts)
h = 45 * 10**(-3)  # m
w_x = 2 * 10**(-3)
w_y = 10 * 10**(-3)
write('build/h.tex', make_SI(h * 1e3, r'\milli\meter', figures=1))
write('build/w_x.tex', make_SI(w_x * 1e3, r'\milli\meter', figures=1))
write('build/w_y.tex', make_SI(w_y * 1e3, r'\milli\meter', figures=1))

Omega = 4 * np.arctan((w_x * w_y) / (2 * h * np.sqrt(4 * h**2 + w_x**2 + w_y**2)))  # Berechne Raumwinkel für eine Pyramide https://de.wikipedia.org/wiki/Raumwinkel
Omega_degree = (Omega * 360 / 2 / np.pi)
write('build/Omega.tex', make_SI(Omega, r'\steradian', figures=2))
write('build/Omega_degree.tex', make_SI(Omega_degree, r'\degree', figures=1))
# x = (7.2 / 11.6) * (1 / (2 * 10**(-6) * 5.91))
# y = (7.2 / 11.6) * (20 * 10**(-6) / 3.03 * 10**14)
F = w_x * w_y
N_T = 19.3 * 6.022 * 10**(23) * 4 * 10**(-5) / 196.96657  # rho[g/cm^3] * N_A *Volumen_Target[cm^3]/Molare Masse [u] https://physik.leech.it/pub/Ex_VI/Folien/SS_13/VL_02.pdf  Folie 4
write('build/F.tex', make_SI(F * 1e6, r'\milli\meter\squared', figures=1))
write('build/N_T.tex', make_SI(N_T, r'', figures=1))

dsigma_domega = Counts * F / (Akt * N_T * Omega)  # https://de.wikipedia.org/wiki/Wirkungsquerschnitt

write('build/Tabelle_t1.tex', make_table([Winkel, Counts, t, dsigma_domega * 1e22], [1, 1, 1, 1]))     # fehlerbehaftete arrays (uarrays) sollten rechts eine 1 bekommen (1 signifikante Stelle)
write('build/Tabelle_t1_texformat.tex', make_full_table(
    caption='Messdaten.',
    label='tab:Winkel',
    source_table='build/Tabelle_t1.tex',
    stacking=[1, 3],              # Hier aufpassen: diese Zahlen bezeichnen diejenigen resultierenden Spaltennummern, die Multicolumns sein sollen
    units=[
        r'$Winkel \:/\: \si{\degree}$',
        r'$Counts$',
        r'$\Delta t \:/\: \si{\second}$',
        r'$\sigma \:/\: 10^{-22} \si{\meter\squared}$'],
    replaceNaN=True,                      # default = false
    replaceNaNby='not a number'))         # default = '-'

dsigma_domega = Counts_norm * F / (Akt * N_T * Omega)  # https://de.wikipedia.org/wiki/Wirkungsquerschnitt
# print(dsigma_domega)


# def Rutherford(w):
#     return (4 * np.pi * epsilon)**(-2) * (z * Z * e / (4 * 5.486 * 10**6))**2 * (1 / (np.sin(np.deg2rad(w))**4))

def Rutherford(w):
    return 1.3 * 10**(-3) * (z * Z / (5.486))**2 * (1 / (np.sin(np.deg2rad(w / 2))**4))  # https://physik.leech.it/pub/Ex_VI/Folien/SS_13/VL_02.pdf Folie 12


w = np.linspace(0, 20, 100)
print(Rutherford(10) * 1e-28)
plt.plot(w, Rutherford(w), 'r-', label='Theorie')
plt.plot(Winkel, dsigma_domega * 1e28, 'bx', label='Messdaten')
plt.xlabel("Winkel / [°] ")
plt.ylabel("Wirkungsquerschnitt / [b/sr]")
plt.ylim(0.5 * min(dsigma_domega * 1e28), 1.3 * max(dsigma_domega * 1e28))
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Wirkungsquerschnitt.pdf')
plt.clf()
##################################################################

vol = np.array([4, 8])
Counts = np.array([1902, 1277])
t = np.array([300, 300])
Counts = unp.uarray(Counts / t, np.sqrt(Counts) / t)
F_d = np.array([2, 4])

N_T = 19.3 * 6.022 * 10**(23) * vol * 10**(-5) / 196.96657  # rho[g/cm^3] * N_A *Volumen_Target[cm^3]/Molare Masse [u] https://physik.leech.it/pub/Ex_VI/Folien/SS_13/VL_02.pdf  Folie 4
dsigma_domega = Counts * F / (Akt * N_T * Omega)  # https://de.wikipedia.org/wiki/Wirkungsquerschnitt

write('build/Tabelle_t2.tex', make_table([F_d, Counts, t, dsigma_domega * 1e22], [1, 1, 1, 1]))     # fehlerbehaftete arrays (uarrays) sollten rechts eine 1 bekommen (1 signifikante Stelle)
write('build/Tabelle_t2_texformat.tex', make_full_table(
    caption='Messung der Counts für zwei verschiedene Foliendicken des gleichen Materials',
    label='tab:divers',
    source_table='build/Tabelle_t2.tex',
    stacking=[1, 3],              # Hier aufpassen: diese Zahlen bezeichnen diejenigen resultierenden Spaltennummern, die Multicolumns sein sollen
    units=[
        r'$\text{Foliendicke} \:/\: \si{\micro\meter}$',
        r'$\text{Counts}$',
        r'$\Delta t \:/\: \si{\second}$',
        r'$\sigma \:/\: 10^{-22} \si{\meter\squared}$'],
    replaceNaN=True,                      # default = false
    replaceNaNby='not a number'))         # default = '-'


##################################################################

N_A = 6.022 * 10**23
Counts_Gold = ufloat(2.85, 0.03)
Counts_Bismut = ufloat(0.35, 0.05)
Counts_Alu = ufloat(0.68, 0.06)

# Counts_Gold = ufloat(614 / 100, 0.05)
# Counts_Bismut = ufloat(603 / 400, 0.05)
# Counts_Alu = ufloat(546 / 400, 0.05)

Z_Gold = 79
Z_Bismut = 83
Z_Alu = 13

d_Gold = 2 * 10**(-6)
d_Bismut = 10**(-6)
d_Alu = 3 * 10**(-6)

rho_Gold = 19.3 * 1000
rho_Bismut = 9.78 * 1000
rho_Alu = 2.7 * 1000

M_Gold = 196.9665 * 1.6e-27
M_Bismut = 208.9804 * 1.6e-27
M_Alu = 26.9815 * 1.6e-27


V_Gold = 10.21 * 10**(-6)  # m^3/mol
V_Bismut = 21.31 * 10**(-6)
V_Alu = 10 * 10**(-6)

F = 20 * 10**(-6)  # m^2

########################################################################################################
n_Gold = rho_Gold * F * d_Gold / M_Gold
n_Bismut = rho_Bismut * F * d_Bismut / M_Bismut
n_Alu = rho_Alu * F * d_Alu / M_Alu

I_Gold = Counts_Gold / (n_Gold * d_Gold)
I_Bismut = Counts_Bismut / (n_Bismut * d_Bismut)
I_Alu = Counts_Alu / (n_Alu * d_Alu)


# def Z(z, a, c):
#     return a * z + c


Z_array = np.array([Z_Gold, Z_Bismut, Z_Alu])
I_norm_array = np.array([I_Gold, I_Bismut, I_Alu])
t_plot = np.linspace(0, 100, 100)
params = ucurve_fit(reg_linear, Z_array, noms(I_norm_array))             # linearer Fit
#params = ucurve_fit(Z, Z_array, noms(I_norm_array))
plt.plot(t_plot, reg_linear(t_plot, *noms(params)), 'r-', label='Fit')

plt.plot(Z_Gold, noms(I_Gold), 'bx', label='Gold')
plt.plot(Z_Bismut, noms(I_Bismut), 'gx', label='Bismut')
plt.plot(Z_Alu, noms(I_Alu), 'mx', label='Alu')
plt.xlabel("Z")
plt.ylabel("normierte Intensität / [1/(sm)]")
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Z.pdf')
plt.clf()

####################################################################################

A_1994 = 330000  # Bq
t = (24 * 12 + 4) * 30 * 24 * 60 * 60  # in Sekunden Oktober 1994 bis Februar 2019
T_halb = 432.2 * 365 * 24 * 60 * 60
lambda_Am = np.log(2) / T_halb
A_2019 = A_1994 * np.exp(-lambda_Am * t)
print(A_2019)

write('build/A_2019.tex', make_SI(A_2019 / 1000, r'\kilo\becquerel', figures=1))

Kugel = 4 * np.pi * 0.101
Blende = 20 * 10**(-6)
Faktor = Blende / Kugel
print(A_2019 * Faktor)

write('build/A_theo.tex', make_SI(A_2019 * Faktor, r'\becquerel', figures=1))

Abweichung_A = (Akt - A_2019 * Faktor) / (A_2019 * Faktor)
write('build/Akt_Abweichung.tex', make_SI(Abweichung_A * 100, r'\percent', figures=1))

Abweichung_d = (d - 2 * 10**(-6)) / (2 * 10**(-6))
write('build/d_Abweichung.tex', make_SI(Abweichung_d * 100, r'\percent', figures=1))

#######################################################################################################

# n_Gold = 20 * 10**(-6) * d_Gold / V_Gold
# n_Bismut = 20 * 10**(-6) * d_Bismut / V_Bismut
# n_Alu = 20 * 10**(-6) * d_Alu / V_Alu


# I_Gold = Counts_Gold / (n_Gold * d_Gold)
# I_Bismut = Counts_Bismut / (n_Bismut * d_Bismut)
# I_Alu = Counts_Alu / (n_Alu * d_Alu)


# # I_Gold = Counts_Gold * V_Gold / (n_Gold* d_Gold)
# # I_Bismut = Counts_Bismut * V_Bismut / (N_A * F * d_Bismut**2)
# # I_Alu = Counts_Alu * V_Alu / (N_A * F * d_Alu**2)

# # I_norm_Gold = I_Gold / (n_Gold * d_Gold)
# # I_norm_Bismut = I_Bismut / (n_Bismut * d_Bismut)
# # I_norm_Alu = I_Alu / (n_Alu * d_Alu)


# def Z(z, a, c):
#     return a * z + c


# Z_array = np.array([Z_Gold, Z_Bismut, Z_Alu])
# I_norm_array = np.array([I_Gold, I_Bismut, I_Alu])
# # t_plot = np.linspace(0, 100, 100)
# #params = ucurve_fit(Z, Z_array, noms(I_norm_array))
# #plt.plot(t_plot, Z(t_plot, *noms(params)), 'r-', label='Fit')

# plt.plot(Z_Gold, noms(I_Gold), 'bx', label='Gold')
# plt.plot(Z_Bismut, noms(I_Bismut), 'gx', label='Bismut')
# plt.plot(Z_Alu, noms(I_Alu), 'mx', label='Alu')
# plt.xlabel("Z")
# plt.ylabel("I / [1/(sm)]")
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/Z.pdf')
# plt.clf()

##############################################################################################################################

# Z = np.array([13, 79, 83])
# n_gold = (19.30 * 1000) / (196.96657 * 1.6e-27)
# n_bismut = (9.78 * 1000) / (208.9 * 1.6e-27)
# n_alu = (2.7 * 1000) / (26.9 * 1.6e-27)
# I_norm = np.array([Counts_Alu / (n_alu * d_Alu), Counts_Gold / (n_gold * d_Gold), Counts_Bismut / (n_bismut * d_Bismut)])


# n = 5.8e22 * 2.5
# X = [13, 79, 83]
# nb = (9.78 * 1000) / (208.9 * 1.6e-27)
# nal = (2.7 * 1000) / (26.9 * 1.6e-27)
# Y = [546 / 400 / (nal * 3e-6), 614 / 100 / n, 603 / 400 / (nb * 1e-6)]


# print(I_norm)
# print(Y)

# plt.plot(Z, noms(I_norm), 'bx', label='Gold')
# plt.xlabel("Z")
# plt.ylabel("I / [1/(sm)]")
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/Z.pdf')
# plt.clf()

################################ FREQUENTLY USED CODE ################################
#
########## IMPORT ##########
# t, U, U_err = np.genfromtxt('messdaten/data.txt', unpack=True)
# t *= 1e-3


########## ERRORS ##########
# R_unc = ufloat(R[0],R[2])
# U = 1e3 * unp.uarray(U, U_err)
# Rx_mean = np.mean(Rx)                 # Mittelwert und syst. Fehler
# Rx_mean_with_error = mean(Rx, 0)      # unp.uarray mit Fehler und Fehler des Mittelwertes, die 0 gibt an, dass in einem R^2 array jeweils die Zeilen gemittelt werden sollen
# Rx_mean_err = MeanError(noms(Rx))     # nur der Fehler des Mittelwertes
#
# Relative Fehler zum späteren Vergleich in der Diskussion
# RelFehler_G = (G_mess - G_lit) / G_lit
# RelFehler_B = (B_mess - B_lit) / B_lit
# write('build/RelFehler_G.tex', make_SI(RelFehler_G*100, r'\percent', figures=1))
# write('build/RelFehler_B.tex', make_SI(RelFehler_B*100, r'\percent', figures=1))


########## CURVE FIT ##########
# def f(t, a, b, c, d):
#     return a * np.sin(b * t + c) + d
#
# params = ucurve_fit(f, t, U, p0=[1, 1e3, 0, 0])   # p0 bezeichnet die Startwerte der zu fittenden Parameter
# params = ucurve_fit(reg_linear, x, y)             # linearer Fit
# params = ucurve_fit(reg_quadratic, x, y)          # quadratischer Fit
# params = ucurve_fit(reg_cubic, x, y)              # kubischer Fit
# a, b = params
# write('build/parameter_a.tex', make_SI(a * 1e-3, r'\kilo\volt', figures=1))       # type in Anz. signifikanter Stellen
# write('build/parameter_b.tex', make_SI(b * 1e-3, r'\kilo\hertz', figures=2))      # type in Anz. signifikanter Stellen

########## MAKE_SI ##########
# make_SI(m_eff_n*1e32, r'\kilo\gramm', figures=1, exp='e-32')
# ufloats, die einen exponenten besitzen, werden als (2.3+/-0.3)e-32 dargestellt, das Problem ist hier der Exponent,
# man sollte sich bei einem Fehler der Form:
# TypeError: non-empty format string passed to object.__format__
# den Wert über print() ausgeben lassen und dann obige Modifikation (links *1e32, und dann exp='e-32') nutzen.
# (hierzu gibt es momentan /Juni2017 ein issue)


########## PLOTTING ##########
# plt.clf()                   # clear actual plot before generating a new one
#
# automatically choosing limits with existing array T1
# t_plot = np.linspace(np.amin(T1), np.amax(T1), 100)
# plt.xlim(t_plot[0]-1/np.size(T1)*(t_plot[-1]-t_plot[0]), t_plot[-1]+1/np.size(T1)*(t_plot[-1]-t_plot[0]))
#
# hard coded limits
# t_plot = np.linspace(-0.5, 2 * np.pi + 0.5, 1000) * 1e-3
#
# standard plotting
# plt.plot(t_plot * 1e3, f(t_plot, *noms(params)) * 1e-3, 'b-', label='Fit')
# plt.plot(t * 1e3, U * 1e3, 'rx', label='Messdaten')
# plt.errorbar(B * 1e3, noms(y) * 1e5, fmt='rx', yerr=stds(y) * 1e5, label='Messdaten')        # mit Fehlerbalken
# plt.xscale('log')                                                                            # logarithmische x-Achse
# plt.xlim(t_plot[0] * 1e3, t_plot[-1] * 1e3)
# plt.xlabel(r'$t \:/\: \si{\milli\second}$')
# plt.ylabel(r'$U \:/\: \si{\kilo\volt}$')
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/aufgabenteil_a_plot.pdf')


########## WRITING TABLES ##########
# IF THERE IS ONLY ONE COLUMN IN A TABLE (workaround):
# a=np.array([Wert_d[0]])
# b=np.array([Rx_mean])
# c=np.array([Rx_mean_err])
# d=np.array([Lx_mean*1e3])
# e=np.array([Lx_mean_err*1e3])
#
# write('build/Tabelle_b.tex', make_table([a,b,c,d,e],[0, 1, 0, 1, 1]))     # fehlerbehaftete arrays (uarrays) sollten rechts eine 1 bekommen (1 signifikante Stelle)
# write('build/Tabelle_b_texformat.tex', make_full_table(
#     caption = 'Messdaten Kapazitätsmessbrücke.',
#     label = 'table:A2',
#     source_table = 'build/Tabelle_b.tex',
#     stacking = [1,2,3,4,5],              # Hier aufpassen: diese Zahlen bezeichnen diejenigen resultierenden Spaltennummern, die Multicolumns sein sollen
#     units = ['Wert',
#     r'$C_2 \:/\: \si{\nano\farad}$',
#     r'$R_2 \:/\: \si{\ohm}$',
#     r'$R_3 / R_4$', '$R_x \:/\: \si{\ohm}$',
#     r'$C_x \:/\: \si{\nano\farad}$'],
#     replaceNaN = True,                      # default = false
#     replaceNaNby = 'not a number'))         # default = '-'
#
# Aufsplitten von Tabellen, falls sie zu lang sind
# t1, t2 = np.array_split(t * 1e3, 2)
# U1, U2 = np.array_split(U * 1e-3, 2)
# write('build/loesung-table.tex', make_table([t1, U1, t2, U2], [3, None, 3, None]))  # type in Nachkommastellen
#
# Verschmelzen von Tabellen (nur Rohdaten, Anzahl der Spalten muss gleich sein)
# write('build/Tabelle_b_composed.tex', make_composed_table(['build/Tabelle_b_teil1.tex','build/Tabelle_b_teil2.tex']))


########## ARRAY FUNCTIONS ##########
# np.arange(2,10)                   # Erzeugt aufwärts zählendes Array von 2 bis 10
# np.zeros(15)                      # Erzeugt Array mit 15 Nullen
# np.ones(15)                       # Erzeugt Array mit 15 Einsen
#
# np.amin(array)                    # Liefert den kleinsten Wert innerhalb eines Arrays
# np.argmin(array)                  # Gibt mir den Index des Minimums eines Arrays zurück
# np.amax(array)                    # Liefert den größten Wert innerhalb eines Arrays
# np.argmax(array)                  # Gibt mir den Index des Maximums eines Arrays zurück
#
# a1,a2 = np.array_split(array, 2)  # Array in zwei Hälften teilen
# np.size(array)                    # Anzahl der Elemente eines Arrays ermitteln
#
# np.delte(A,3)                     # liefert ein Array, in dem der Eintrag mit Index 3 des arrays
# A gelöscht wurde (nachfolgende indices werden aufgerückt)


########## ARRAY INDEXING ##########
# y[n - 1::n]                       # liefert aus einem Array jeden n-ten Wert als Array


########## DIFFERENT STUFF ##########
# R = const.physical_constants["molar gas constant"]      # Array of value, unit, error, siehe https://docs.scipy.org/doc/scipy-0.19.0/reference/constants.html
# search_replace_within_file('build/Tabelle_test.tex','find me','found you')    # Selbsterklärend


# convenience file writing for standard make Files
f = open('build/.pysuccess', 'w')
f.write('MarktkraM')
f.close()
