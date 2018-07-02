##################################################### Import system libraries ######################################################
import matplotlib as mpl
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
from scipy.optimize import curve_fit
##############Arrays die nachher benötigt werden########
t_0=np.array([])
U_1=np.array([])
U_2=np.array([])
U_fit=np.array([])

##################################Berechnung von T1#########################################
#Laden der Daten
tau, A = np.genfromtxt("messdaten/t1_bestimmung.txt", unpack=True)

def fitfunction(tau, a, T1, b):
     return -a * ( 1 -( 2 * np.exp(-tau/T1)))+b

params_1, covariance_1 = curve_fit( fitfunction, tau, A , maxfev=900000000)
errors_1 = np.sqrt(np.diag(covariance_1))
#Plotten der Daten
x_plot = np.linspace(0, 5000 ,10000)
plt.plot(x_plot, fitfunction(x_plot, *params_1), "r-", label="Fit")
plt.plot( tau, A, ".c",label="Messwerte zur Bestimmung von Tau_1")
plt.xlabel('tau in ms')
plt.ylabel('A in mV')
plt.legend(loc="best")
plt.grid()
plt.savefig("build/fit_tau_1.pdf")
plt.clf()



###############################################Tau2 Bestimmen
t_xaxsis, U1, U2 =np.genfromtxt('messdaten/scope_t2.csv', delimiter=',',comments="#", unpack=True)
U2=-1*U2
def fitfunction_tau2(tau_values, a, b, c):
     return -a * ( np.exp(-tau_values/b))+c

params_2, covariance_2 = curve_fit( fitfunction_tau2, t_xaxsis, U2 , maxfev=900000000)
errors_2 = np.sqrt(np.diag(covariance_2))
plt.clf()

plt.plot(t_xaxsis,U2,".k")
x_plot_tau2=np.linspace(0,1.905,10000)
plt.plot(x_plot_tau2, fitfunction_tau2(x_plot_tau2, *params_2), "r-", label="Fit_Tau2")
plt.xlabel('tau in ms')
plt.ylabel('A in V')
plt.grid()
plt.savefig("build/U_fit_Tau2.pdf")
plt.clf()



###########Bestimmung der Diffusionskonstante D #################
#es sieht so aus als ist irgenwas falsch gepolt, da bei uns vorzeichen und ähnliches Vertauscht sind

td, Udiff = np.genfromtxt("messdaten/diff_const.txt", unpack=True)
def fitfunction_D(td,d, a,c,b):
     return -a * (np.exp(-td/d)*np.exp(-b* td**3))+c
Udiff=Udiff
p0=[-10,700, 10,0.0002]
params_3, covariance_3 = curve_fit( fitfunction_D, td, Udiff
,p0=p0
, maxfev=999000000)
errors_3 = np.sqrt(np.diag(covariance_2))

x_range=np.linspace(5,31,10000)
plt.plot(td,Udiff,".r")
plt.plot(x_range, fitfunction_D(x_range, *params_3), "r-", label="Fit_D")
plt.plot( td, Udiff, ".b",label="Messwerte zur Bestimmung von D")
plt.xlabel('tau in ms')
plt.grid()
plt.ylabel('A in mV')
plt.legend(loc="best")
plt.savefig("build/plot_diff.pdf")
plt.clf()


##Bestimmung der Halbwertsbreite
t_w, U1, U2 =np.genfromtxt('messdaten/scope_1.csv', delimiter=',', unpack=True)
plt.plot(t_w, U2, ".c")
plt.grid()
plt.axvline(x=0.020127177, c="black")
plt.axvline(x=0.020010, c="black")
plt.axhline(y=0.56177, c="blue")
plt.axhline(y=0.28085, c="green")
plt.axhline(y=0.22085, xmin=0.020, xmax=0.0201,c="red")
plt.axhline(y=0.0, c="blue")
plt.savefig("build/halbw.pdf")
plt.clf()



####Plot der Parameter um diese abzuschätzen
#x_range_plot=np.linspace(0,31,1000)
#a,c,b,d=p0
#values=a * (np.exp(-x_range_plot/T2)*np.exp(-b* x_range_plot**3))+c
#plt.plot(x_range_plot, values, ".k")
#plt.savefig("build/paraschäter.pdf")
#plt.clf()




#Print der Parameter T1
print("Parameter der Fits, für die T1 Bestimmung")
print("a=0", params_1[0], "±", errors_1[0])
print("T_1=0", params_1[1], "±", errors_1[1])

#Print der Parameter T2
print("Parameter der Fits zur Bestimmung von T2")
print("U_0=0", params_2[0], "±", errors_2[0])
print("T_2=0", params_2[1], "±", errors_2[1])
print("c", params_2[2], "±", errors_2[2])

#Print der Parameter  D
print("Parameter der Fits von D")
print("a", params_3[0], "±", errors_3[0])
print("b", params_3[1], "±", errors_3[1])
print("c", params_3[2], "±", errors_3[2])






f = open('build/.pysuccess', 'w')
f.write('MarktkraM')
f.close()
