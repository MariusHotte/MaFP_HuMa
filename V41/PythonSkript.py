##################################################### Import system libraries ######################################################
import matplotlib
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)
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
import os, sys, inspect
# realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"python_custom_scripts")))
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

Metall_cm =np.array([4,5.7,7.1,8.4,9.6,10.9,12.2,13.8,16.4]) # Diese Angaben sind in cm. 180Grad entsprechen auf dem Streifen 18cm Metall_cm =np.array([4,5.7,7.1,8.4,9.6,10.9,12.2,13.8,16.4])
Salz_cm =np.array([2.75,3.9,4.8,6.3,7.0,8.3,8.9,9.55,10.2,11.55,12.2,14.9,16.5])
m_bcc= np.array([1.41,2.00,2.45,2.83,3.16,3.46,3.74,4.00])#Metall m h k l bcc
lam=1.54093*10**(-10) # Wellenlänge des Röntgenstrahls in meter
Miller_diamant=np.array([111,220,311,400,331,422,333,440,531]) ##Alles 4 fürs metall
Miller_bcc=np.array([110,200,211,220,310,222,321,400,330])
Miller_fcc=np.array([111,200,220,311,222,400,331,420,422])
Miller_sc=np.array([100,110,111,200,210,211,220,221,310])

m_diamant=np.array([])
for i in range(len(Miller_diamant)):
	h=int(list(str(Miller_diamant[i]))[0])
	k=int(list(str(Miller_diamant[i]))[1])
	l=int(list(str(Miller_diamant[i]))[2])
	x=np.sqrt((h**2)+(k**2)+(l**2))
	m_diamant=np.append(m_diamant,[x])
m_diamant_norm=m_diamant/m_diamant[0]

m_bcc=np.array([])
for i in range(len(Miller_bcc)):
	h=int(list(str(Miller_bcc[i]))[0])
	k=int(list(str(Miller_bcc[i]))[1])
	l=int(list(str(Miller_bcc[i]))[2])
	x=np.sqrt((h**2)+(k**2)+(l**2))
	m_bcc=np.append(m_bcc,[x])
m_bcc_norm=m_bcc/m_bcc[0]

m_fcc=np.array([])
for i in range(len(Miller_fcc)):
	h=int(list(str(Miller_fcc[i]))[0])
	k=int(list(str(Miller_fcc[i]))[1])
	l=int(list(str(Miller_fcc[i]))[2])
	x=np.sqrt((h**2)+(k**2)+(l**2))
	m_fcc=np.append(m_fcc,[x])
m_fcc_norm=m_fcc/m_fcc[0]


m_sc=np.array([])
for i in range(len(Miller_sc)):
	h=int(list(str(Miller_sc[i]))[0])
	k=int(list(str(Miller_sc[i]))[1])
	l=int(list(str(Miller_sc[i]))[2])
	x=np.sqrt((h**2)+(k**2)+(l**2))
	m_sc=np.append(m_sc,[x])
m_sc_norm=m_sc/m_sc[0]



###Metall##################
Metall_cm = unp.uarray(Metall_cm, 0.1)   # nehme einen Ablesefehler von einem Millimeter an
theta_strich= Metall_cm/5.73
theta=theta_strich/2

d=lam/(2*unp.sin(theta))
d_verhaeltnis=(d[0]/d)
# print(d_verhaeltnis)

write('build/Tabelle_Metall_Rohdaten.tex', make_table([Metall_cm,theta,d*10**10],[1,1,1,1,1,1,]))     # fehlerbehaftete arrays (uarrays) sollten rechts eine 1 bekommen (1 signifikante Stelle)
write('build/Tabelle_Metall_Rohdaten_texformat.tex', make_full_table(
    caption = 'Messdaten der Metallprobe.',
    label = 'table:A1',
    source_table = 'build/Tabelle_Metall_Rohdaten.tex',
    stacking = [0,1,2],              # Hier aufpassen: diese Zahlen bezeichnen diejenigen resultierenden Spaltennummern, die Multicolumns sein sollen
    units = [
    r'$s\:/\: \si{\centi\meter}$',
    r'$\theta\:/\: \si{\centi\meter}$',
    r'$d \:/\: \si{\angstrom}$'],
    replaceNaN = True,                      # de\fault = false
    replaceNaNby = 'not a number'))         # default = '-'


write('build/Tabelle_Metall_Kristall.tex', make_table([Miller_sc,m_sc_norm,Miller_bcc,m_bcc_norm,Miller_fcc,m_fcc_norm,Miller_diamant,m_diamant_norm,d_verhaeltnis],[0,2,0,2,0,2,0,2,1]))     # fehlerbehaftete arrays (uarrays) sollten rechts eine 1 bekommen (1 signifikante Stelle)
write('build/Tabelle_Metall_Kristall_texformat.tex', make_full_table(
    caption = 'Messdaten der Metallprobe.',
    label = 'table:A2',
    source_table = 'build/Tabelle_Metall_Kristall.tex',
    stacking = [8],              # Hier aufpassen: diese Zahlen bezeichnen diejenigen resultierenden Spaltennummern, die Multicolumns sein sollen
    units = [
    r'$\text{sc}$',
    r'$\left(\frac{m_i}{m_1}\right)_{sc}$',
    r'$\text{bcc}$',
    r'$\left(\frac{m_i}{m_1}\right)_{bcc}$',
    r'$\text{fcc}$',
    r'$\left(\frac{m_i}{m_1}\right)_{fcc}$',
    r'$\text{diamant}$',
    r'$\left(\frac{m_i}{m_1}\right)_{dia}$',
    r'$\frac{d_1}{d_i}$'],
    replaceNaN = True,                      # de\fault = false
    replaceNaNby = 'not a number'))         # default = '-'

a=noms(d*m_bcc)
x=noms((unp.cos(theta))**2)
params = ucurve_fit(reg_linear, x, a*10**10)             # linearer Fit
t_plot = np.linspace(0, 1, 2)
plt.plot(t_plot, reg_linear(t_plot, *noms(params)), 'b-', label='Fit')
plt.plot(x, a*10**10, 'rx', label='Messdaten')
#plt.xlabel(r'$a \:/\: \si{\angstrom}$')
#plt.ylabel(r'$\cos^2{(theta)}$')
plt.legend(loc='best')
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Metall.pdf')
plt.clf()

write('build/a_Metall.tex', make_SI(params[1], r'\angstrom', figures=2))
write('build/b_Metall.tex', make_SI(params[0], r'\angstrom', figures=2))

lit_wert=3.30
prozent=1-params[0]/lit_wert
write('build/prozent_metall.tex', make_SI(prozent, r'\percent', figures=1))
# ###Salz##################

theta_strich= Salz_cm/5.73
theta=theta_strich/2

d=lam/(2*np.sin(theta))
# d_verhaeltnis=(d[0]/d)
# print(d_verhaeltnis)
import cmath
ii= complex(0,1)

print('Jetzt geht es los!')
Miller_zinkblende=np.array([])
Miller=([100,110,111,200,210,211,220,221,222,300,310,311,320,321,322,330,331,332,333,400,410,411,420,421,422,430,431,432,433,440,441,442,443,444])
for i in range(len(Miller)):
	h=int(list(str(Miller[i]))[0])
	k=int(list(str(Miller[i]))[1])
	l=int(list(str(Miller[i]))[2])
	atom_a= cmath.exp(-2*np.pi*ii*(0*h+0*k+0*l)) + cmath.exp(-2*np.pi*ii*(0.5*h+0.5*k+0*l)) + cmath.exp(-2*np.pi*ii*(0*h+0.5*k+0.5*l)) +cmath.exp(-2*np.pi*ii*(0.5*h+0*k+0.5*l))
	atom_b= cmath.exp(-2*np.pi*ii*(0.25*h+0.25*k+0.25*l)) + cmath.exp(-2*np.pi*ii*(0.75*h+0.75*k+0.25*l)) + cmath.exp(-2*np.pi*ii*(0.25*h+0.75*k+0.75*l)) +cmath.exp(-2*np.pi*ii*(0.75*h+0.25*k+0.75*l))
	if abs(atom_a.real)>10**(-3) or abs(atom_a.imag)>10**(-3) or abs(atom_b.real)>10**(-3) or abs(atom_b.imag)>10**(-3):
		Miller_zinkblende=np.append(Miller_zinkblende,Miller[i])
m_zinkblende=np.array([])
for i in range(len(Miller_zinkblende)):
    h=int(list(str(Miller_zinkblende[i]))[0])
    k=int(list(str(Miller_zinkblende[i]))[1])
    l=int(list(str(Miller_zinkblende[i]))[2])
    x=np.sqrt((h**2)+(k**2)+(l**2))
    m_zinkblende=np.append(m_zinkblende,[x])


Miller_steinsalz=np.array([])
for i in range(len(Miller)):
	h=int(list(str(Miller[i]))[0])
	k=int(list(str(Miller[i]))[1])
	l=int(list(str(Miller[i]))[2])
	atom_a= cmath.exp(-2*np.pi*ii*(0*h+0*k+0*l)) + cmath.exp(-2*np.pi*ii*(0.5*h+0.5*k+0*l)) + cmath.exp(-2*np.pi*ii*(0*h+0.5*k+0.5*l)) +cmath.exp(-2*np.pi*ii*(0.5*h+0*k+0.5*l))
	atom_b= cmath.exp(-2*np.pi*ii*(0.5*h+0.5*k+0.5*l)) + cmath.exp(-2*np.pi*ii*(1*h+1*k+0.5*l)) + cmath.exp(-2*np.pi*ii*(1*h+0.5*k+1*l)) +cmath.exp(-2*np.pi*ii*(0.5*h+1*k+1*l))
	if abs(atom_a.real)>10**(-3) or abs(atom_a.imag)>10**(-3) or abs(atom_b.real)>10**(-3) or abs(atom_b.imag)>10**(-3):
		Miller_steinsalz=np.append(Miller_steinsalz,Miller[i])
m_steinsalz=np.array([])
for i in range(len(Miller_steinsalz)):
    h=int(list(str(Miller_steinsalz[i]))[0])
    k=int(list(str(Miller_steinsalz[i]))[1])
    l=int(list(str(Miller_steinsalz[i]))[2])
    x=np.sqrt((h**2)+(k**2)+(l**2))
    m_steinsalz=np.append(m_steinsalz,[x])

Miller_chlorid=np.array([])
for i in range(len(Miller)):
	h=int(list(str(Miller[i]))[0])
	k=int(list(str(Miller[i]))[1])
	l=int(list(str(Miller[i]))[2])
	atom_a= cmath.exp(-2*np.pi*ii*(0*h+0*k+0*l))+ cmath.exp(-2*np.pi*ii*(0.5*h+0.5*k+0.5*l))
	if abs(atom_a.real)>10**(-3) or abs(atom_a.imag)>10**(-3):
		Miller_chlorid=np.append(Miller_chlorid,Miller[i])
m_chlorid=np.array([])
for i in range(len(Miller_chlorid)):
    h=int(list(str(Miller_chlorid[i]))[0])
    k=int(list(str(Miller_chlorid[i]))[1])
    l=int(list(str(Miller_chlorid[i]))[2])
    x=np.sqrt((h**2)+(k**2)+(l**2))
    m_chlorid=np.append(m_chlorid,[x])


Miller_caesium=np.array([100,110,111,200,221,210,211,220,300,310,311,222,320,321,330,400,410,411,331,420,421,332,422,333,440,442,500,510,324,521])
m_caesium=np.array([])
for i in range(len(Miller_caesium)):
	h=int(list(str(Miller_caesium[i]))[0])
	k=int(list(str(Miller_caesium[i]))[1])
	l=int(list(str(Miller_caesium[i]))[2])
	x=np.sqrt((h**2)+(k**2)+(l**2))
	m_caesium=np.append(m_caesium,[x])

print(m_caesium,m_chlorid, m_zinkblende, m_steinsalz)
# a_0=(lam/2*np.sin(theta[0]))*m_caesium

# n = 13
# m = 30
# a = [0] * n
# for i in range(n):
#     a[i] = [0] * m
# for j in range(len(theta)):
# 	for i in range(len(m_caesium)):
# 		a[j][i]=(lam/2*np.sin(theta[j]))*m_caesium[i]*10**10
# b_1=np.array(a[0])
# b_2=np.array(a[1])

# print('Hier')
# print(len(b_1))
# print(len(Miller_caesium))



# write('build/Tabelle_1.tex', make_table([Miller_caesium,b_1,b_2],[1,1,1]))    # fehlerbehaftete arrays (uarrays) sollten rechts eine 1 bekommen (1 signifikante Stelle)
# write('build/Tabelle_1_texformat.tex', make_full_table(
#     caption = 'Chaesiumchlorid.',
#     label = 'table:A3',
#     source_table = 'build/Tabelle_1.tex',
#     stacking = [],              # Hier aufpassen: diese Zahlen bezeichnen diejenigen resultierenden Spaltennummern, die Multicolumns sein sollen
#     units = [
#     r'$a_1$',
#     r'$a_1$',
#     r'$a_2$',],
#     replaceNaN = True,                      # default = false
#     replaceNaNby = 'not a number'))         # default = '-'







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
## Relative Fehler zum späteren Vergleich in der Diskussion
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
## automatically choosing limits with existing array T1
# t_plot = np.linspace(np.amin(T1), np.amax(T1), 100)
# plt.xlim(t_plot[0]-1/np.size(T1)*(t_plot[-1]-t_plot[0]), t_plot[-1]+1/np.size(T1)*(t_plot[-1]-t_plot[0]))
#
## hard coded limits
# t_plot = np.linspace(-0.5, 2 * np.pi + 0.5, 1000) * 1e-3
#
## standard plotting
# plt.plot(t_plot * 1e3, f(t_plot, *noms(params)) * 1e-3, 'b-', label='Fit')
# plt.plot(t * 1e3, U * 1e3, 'rx', label='Messdaten')
## plt.errorbar(B * 1e3, noms(y) * 1e5, fmt='rx', yerr=stds(y) * 1e5, label='Messdaten')        # mit Fehlerbalken
## plt.xscale('log')                                                                            # logarithmische x-Achse
# plt.xlim(t_plot[0] * 1e3, t_plot[-1] * 1e3)
# plt.xlabel(r'$t \:/\: \si{\milli\second}$')
# plt.ylabel(r'$U \:/\: \si{\kilo\volt}$')
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/aufgabenteil_a_plot.pdf')


########## WRITING TABLES ##########
### IF THERE IS ONLY ONE COLUMN IN A TABLE (workaround):
## a=np.array([Wert_d[0]])
## b=np.array([Rx_mean])
## c=np.array([Rx_mean_err])
## d=np.array([Lx_mean*1e3])
## e=np.array([Lx_mean_err*1e3])
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
## Aufsplitten von Tabellen, falls sie zu lang sind
# t1, t2 = np.array_split(t * 1e3, 2)
# U1, U2 = np.array_split(U * 1e-3, 2)
# write('build/loesung-table.tex', make_table([t1, U1, t2, U2], [3, None, 3, None]))  # type in Nachkommastellen
#
## Verschmelzen von Tabellen (nur Rohdaten, Anzahl der Spalten muss gleich sein)
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


## convenience file writing for standard make Files
f = open('build/.pysuccess', 'w')
f.write('MarktkraM')
f.close()
