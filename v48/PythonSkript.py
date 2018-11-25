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

# u = np.array([0.249, 0.198, 0.37, 0.174, 0.134])
# z = np.array([0.05341, 0.0643, 0.05284, 0.0644, 0.04912])
# for i in range(0, len(z)):
#     U = ufloat(u[i], 0.02 * u[i])
#     I = ufloat(z[i], 0.02 * z[i])
#     print(U / I)
################################ FREQUENTLY USED CODE ################################
# T  A  Offset I=0.3*10^(-9)
Temp, Strom = np.genfromtxt('messdaten/erste_messung_mod.txt', unpack=True)
Temp = Temp + 273.15
Offset = 0.3 * 10**(-11)
Strom = Strom * (-1) * 10**(-11)
Temp_err = np.zeros(len(Temp))
Temp_all = unp.uarray(Temp, Temp_err)
Strom_err = np.zeros(len(Temp))
Strom_all = unp.uarray(Strom, Strom_err)

t_plot = np.linspace(0, len(Temp), len(Temp))
background_Temp = Temp[15:20]
background_Temp = np.append(background_Temp, Temp[35:41])
background_Strom = Strom[15:20]
background_Strom = np.append(background_Strom, Strom[35:41])

plt.plot(t_plot, Temp, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_Time.pdf')
plt.clf()


plt.plot(Temp, Strom, 'rx', label='Messdaten')
#plt.plot(background_Temp, background_Strom, 'bx', label='Background')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current.pdf')
plt.clf()


def background(T, a, b, c):
    return a * np.exp(b * T) + c


params = ucurve_fit(background, background_Temp, background_Strom, p0=[13 * 10**(-15), 0.03, 225])   # p0 bezeichnet die Startwerte der zu fittenden Parameter

write('build/1_peak1_untergrundparams_a.tex', make_SI(params[0] * 10**12, r'\pico\ampere', figures=2))
write('build/1_peak1_untergrundparams_b.tex', make_SI(params[1], r'\kelvin', figures=2))
write('build/1_peak1_untergrundparams_c.tex', make_SI(params[2] * 10**9, r'\nano\ampere', figures=2))
#######################################################################################################

# korrigiere ersten Peak


Temp_peak1 = Temp[20:35]
Strom_peak1 = Strom_all[20:35]
Temp_peak2 = Temp[41:len(Temp) - 3]
Strom_peak2 = Strom_all[41:len(Temp) - 3]


plt.plot(Temp, Strom, 'rx', label='Messdaten')
plt.plot(background_Temp, np.abs(background_Strom), 'yx', label='Untergrund')
plt.plot(Temp_peak1, noms(Strom_peak1), 'gx', label='Peak1')
plt.plot(Temp_peak2, noms(Strom_peak2), 'mx', label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_background_peak_withoutfit.pdf')
plt.clf()


t_plot = np.linspace(213.15, 333.15, 1000)
plt.plot(t_plot, background(t_plot, *noms(params)), 'b-', label='Fit')
plt.plot(Temp, Strom, 'rx', label='Messdaten')
plt.plot(background_Temp, np.abs(background_Strom), 'yx', label='Untergrund')
plt.plot(Temp_peak1, noms(Strom_peak1), 'gx', label='Peak1')
plt.plot(Temp_peak2, noms(Strom_peak2), 'mx', label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_background_peak.pdf')
plt.clf()


def background(T, a, b, c):
    return a * unp.exp(b * T) + c


for i in range(len(Temp_peak1)):
    Strom_peak1[i] = Strom_peak1[i] - background(Temp_peak1[i], *params)

for i in range(len(Temp_peak2)):
    Strom_peak2[i] = Strom_peak2[i] - background(Temp_peak2[i], *params)

plt.errorbar(Temp_peak1, noms(Strom_peak1), fmt='rx', yerr=stds(Strom_peak1), label='Peak 1')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak.pdf')
plt.clf()

plt.errorbar(Temp_peak2, noms(Strom_peak2), fmt='rx', yerr=stds(Strom_peak2), label='Peak 2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak2.pdf')
plt.clf()


Strom_peak1_log = unp.log(Strom_peak1)
plt.errorbar(1 / (Temp_peak1[0:10]), noms(Strom_peak1_log[0:10]), fmt='rx', yerr=stds(Strom_peak1_log[0:10]), label='Peak')
# plt.plot(1 / (Temp_peak1[0:10]), np.log(Strom_peak1[0:10]), 'gx', label='Peak')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak_log.pdf')
plt.clf()

Strom_peak2_log = unp.log(Strom_peak2)
plt.errorbar(1 / (Temp_peak2[0:15]), noms(Strom_peak2_log[0:15]), fmt='rx', yerr=stds(Strom_peak2_log[0:15]), label='Peak')
#plt.plot(1 / (Temp_peak2[0:9]), np.log(Strom_peak2[0:9]), 'gx', label='Peak')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak2_log.pdf')
plt.clf()

k_B = 8.6173 * 10**(-5)


def gerade(t, a, b):
    return a * t + b


time = np.linspace(0, len(Temp_peak1), len(Temp_peak1))
params_temp_peak1 = ucurve_fit(gerade, time, Temp_peak1)
t_plot = np.linspace(0, 16, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_temp_peak1)), 'b-', label='Fit')
plt.plot(time, Temp_peak1, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_Time_Peak1.pdf')
plt.clf()

write('build/1_Heizrate_peak1.tex', make_SI(params_temp_peak1[0], r'\kelvin\per\minute', figures=1))

time = np.linspace(0, len(Temp_peak2), len(Temp_peak2))
params_temp_peak2 = ucurve_fit(gerade, time, Temp_peak2)
t_plot = np.linspace(0, 21, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_temp_peak2)), 'b-', label='Fit')
plt.plot(time, Temp_peak2, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_Time_Peak2.pdf')
plt.clf()

write('build/1_Heizrate_peak2.tex', make_SI(params_temp_peak2[0], r'\kelvin\per\minute', figures=1))


def gerade(invT, a, b):
    return (-a / k_B) * invT + b


params = ucurve_fit(gerade, 1 / (Temp_peak1[0:10]), noms(Strom_peak1_log[0:10]))   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(0.0038, 0.0041, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params)), 'b-', label='Fit')
plt.errorbar(1 / (Temp_peak1[0:10]), noms(Strom_peak1_log[0:10]), fmt='rx', yerr=stds(Strom_peak1_log[0:10]), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak_log_fit.pdf')
plt.clf()

tau_0 = (1 / params_temp_peak1[0]) * unp.exp(-params[0] / (k_B * Temp_peak1[10])) * k_B * Temp_peak1[10]**2 / params[0]
write('build/1_tau_peak1_variante1.tex', make_SI(tau_0 * 10**15, r'\pico\second', figures=2))       # type in Anz. signifikanter Stellen
write('build/1_tau_peak1_variante1_temp.tex', make_SI(Temp_peak1[10], r'\kelvin', figures=2))       # type in Anz. signifikanter Stellen

write('build/1_Energie_peak1_variante1.tex', make_SI(params[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen


params = ucurve_fit(gerade, 1 / (Temp_peak2[0:15]), noms(Strom_peak2_log[0:15]))   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(0.003173, 0.0035, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params)), 'b-', label='Fit')
plt.errorbar(1 / (Temp_peak2[0:15]), noms(Strom_peak2_log[0:15]), fmt='rx', yerr=stds(Strom_peak2_log[0:15]), label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak2_log_fit.pdf')
plt.clf()

tau_0 = (1 / params_temp_peak2[0]) * unp.exp(-params[0] / (k_B * Temp_peak2[15])) * k_B * Temp_peak2[15]**2 / params[0]
write('build/1_tau_peak2_variante1.tex', make_SI(tau_0 * 10**15, r'\pico\second', figures=2))       # type in Anz. signifikanter Stellen
write('build/1_tau_peak2_variante1_temp.tex', make_SI(Temp_peak2[15], r'\kelvin', figures=2))       # type in Anz. signifikanter Stellen

write('build/1_Energie_peak2_variante1.tex', make_SI(params[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen
a = params[0]

# plt.plot(Temp_peak1, Strom_peak1, 'gx', label='Peak')
# plt.xlabel(r'$Temp \:/\: \: [K]$')
# plt.ylabel(r'$I \:/\: [A]$')
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/help.pdf')
# plt.clf()


#Strom_peak1_integral = np.ones(len(Temp_peak1))
Strom_peak1_integral = unp.uarray(np.zeros(len(Temp_peak1)), np.zeros(len(Temp_peak1)))
Strom_peak2_integral = unp.uarray(np.zeros(len(Temp_peak2)), np.zeros(len(Temp_peak2)))

for j in range(1, len(Temp_peak1)):
    help = ufloat(0, 0)
    for i in range(j, len(Temp_peak1)):
        help = help + (Temp_peak1[i] - Temp_peak1[i - 1]) * (Strom_peak1[i] + Strom_peak1[i - 1]) / 2   # Spektralverfahren
    Strom_peak1_integral[j - 1] = help

for j in range(1, len(Temp_peak2)):
    help = ufloat(0, 0)
    for i in range(j, len(Temp_peak2)):
        help = help + (Temp_peak2[i] - Temp_peak2[i - 1]) * (Strom_peak2[i] + Strom_peak2[i - 1]) / 2   # Trapezverfahren
    Strom_peak2_integral[j - 1] = help

ln_Strom_peak1_mod = unp.log(Strom_peak1_integral[0:len(Strom_peak1_integral) - 1]) - unp.log(Strom_peak1[0:len(Strom_peak1) - 1])
ln_Strom_peak2_mod = unp.log(Strom_peak2_integral[0:len(Strom_peak2_integral) - 1]) - unp.log(Strom_peak2[0:len(Strom_peak2) - 1])


#plt.plot(Temp_peak1[0:15], ln_Strom_peak1_mod[0:15], 'bx', label='Peak')
plt.errorbar(Temp_peak1[0:len(Temp_peak1) - 1], noms(ln_Strom_peak1_mod), fmt='rx', yerr=stds(ln_Strom_peak1_mod), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$ln(I) \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak_log_2.pdf')
plt.clf()

plt.errorbar(Temp_peak2[0:len(Temp_peak2) - 1], noms(ln_Strom_peak2_mod), fmt='rx', yerr=stds(ln_Strom_peak2_mod), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$ln(I) \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak2_log_2.pdf')
plt.clf()


def gerade(T, a, b):
    return (a / (T * k_B)) + b


params_peak1 = ucurve_fit(gerade, Temp_peak1[0:len(Temp_peak1) - 1], ln_Strom_peak1_mod)   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(243, 273, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_peak1)), 'b-', label='Fit')
plt.errorbar(Temp_peak1[0:len(Temp_peak1) - 1], noms(ln_Strom_peak1_mod), fmt='rx', yerr=stds(ln_Strom_peak1_mod), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$\Omega$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak_log_3.pdf')
plt.clf()

params_peak2 = ucurve_fit(gerade, Temp_peak2[0:len(Temp_peak2) - 1], ln_Strom_peak2_mod)   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(281, 325, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_peak2)), 'b-', label='Fit')
plt.errorbar(Temp_peak2[0:len(Temp_peak2) - 1], noms(ln_Strom_peak2_mod), fmt='rx', yerr=stds(ln_Strom_peak2_mod), label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$\Omega$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/1_Temp_current_peak2_log_3.pdf')
plt.clf()

t_0 = (unp.exp(params_peak1[1]) / params_temp_peak1[0])
write('build/1_tau_peak1_variante2.tex', make_SI(t_0 * 10**18, r'\atto\second', figures=2))       # type in Anz. signifikanter Stellen
write('build/1_Energie_peak1_variante2.tex', make_SI(params_peak1[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen


t_0 = (unp.exp(params_peak2[1]) / params_temp_peak2[0])
write('build/1_tau_peak2_variante2.tex', make_SI(t_0 * 10**18, r'\atto\second', figures=2))       # type in Anz. signifikanter Stellen
write('build/1_Energie_peak2_variante2.tex', make_SI(params_peak2[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen
b = params_peak2[0]
write('build/unterschied.tex', make_SI(1 - a / b, r'\percent', figures=2))       # type in Anz. signifikanter Stellen
#####################################################################################################################################################################
#####################################################################################################################################################################
Temp, Strom = np.genfromtxt('messdaten/zweite.txt', unpack=True)
Temp = Temp + 273.15
Offset = 0.3 * 10**(-11)
Strom = Strom * 10**(-11)
Temp_err = np.zeros(len(Temp))
Temp_all = unp.uarray(Temp, Temp_err)
Strom_err = np.zeros(len(Temp))
Strom_all = unp.uarray(Strom, Strom_err)

t_plot = np.linspace(0, len(Temp), len(Temp))
background_Temp = Temp[:20]
background_Temp = np.append(background_Temp, Temp[42:47])
background_Strom = Strom[:20]
background_Strom = np.append(background_Strom, Strom[42:47])

plt.plot(t_plot, Temp, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_Time.pdf')
plt.clf()


plt.plot(Temp, Strom, 'rx', label='Messdaten')
plt.plot(background_Temp, background_Strom, 'bx', label='Background')
plt.plot(Temp[20:42], Strom[20:42], 'gx', label='Peak1')
plt.plot(Temp[55:], Strom[55:], 'yx', label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current.pdf')
plt.clf()

##################################################################


def background(T, a, b, c):
    return a * np.exp(b * T) + c


params = ucurve_fit(background, background_Temp, background_Strom, p0=[13 * 10**(-15), 0.03, 225])   # p0 bezeichnet die Startwerte der zu fittenden Parameter

write('build/2_peak1_untergrundparams_a.tex', make_SI(params[0] * 10**12, r'\pico\ampere', figures=2))
write('build/2_peak1_untergrundparams_b.tex', make_SI(params[1], r'\kelvin', figures=2))
write('build/2_peak1_untergrundparams_c.tex', make_SI(params[2] * 10**9, r'\nano\ampere', figures=2))

Temp_peak1 = Temp[20:42]
Strom_peak1 = Strom_all[20:42]
Temp_peak2 = Temp[47:len(Temp) - 3]
Strom_peak2 = Strom_all[47:len(Temp) - 3]


t_plot = np.linspace(213.15, 333.15, 1000)
plt.plot(t_plot, background(t_plot, *noms(params)), 'b-', label='Fit')
plt.plot(Temp, Strom, 'rx', label='Messdaten')
plt.plot(background_Temp, np.abs(background_Strom), 'yx', label='Untergrund')
plt.plot(Temp_peak1, noms(Strom_peak1), 'gx', label='Peak1')
plt.plot(Temp_peak2, noms(Strom_peak2), 'mx', label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_background_peak.pdf')
plt.clf()


def background(T, a, b, c):
    return a * unp.exp(b * T) + c


for i in range(len(Temp_peak1)):
    Strom_peak1[i] = Strom_peak1[i] - background(Temp_peak1[i], *params)

for i in range(len(Temp_peak2)):
    Strom_peak2[i] = Strom_peak2[i] - background(Temp_peak2[i], *params)

plt.errorbar(Temp_peak1, noms(Strom_peak1), fmt='rx', yerr=stds(Strom_peak1), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak.pdf')
plt.clf()

plt.errorbar(Temp_peak2, noms(Strom_peak2), fmt='rx', yerr=stds(Strom_peak2), label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$I \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak2.pdf')
plt.clf()


Strom_peak1_log = unp.log(Strom_peak1)
plt.errorbar(1 / (Temp_peak1[0:12]), noms(Strom_peak1_log[0:12]), fmt='rx', yerr=stds(Strom_peak1_log[0:12]), label='Peak')
# plt.plot(1 / (Temp_peak1[0:10]), np.log(Strom_peak1[0:10]), 'gx', label='Peak')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak_log.pdf')
plt.clf()

Strom_peak2_log = unp.log(Strom_peak2)
plt.errorbar(1 / (Temp_peak2[0:22]), noms(Strom_peak2_log[0:22]), fmt='rx', yerr=stds(Strom_peak2_log[0:22]), label='Peak')
#plt.plot(1 / (Temp_peak2[0:9]), np.log(Strom_peak2[0:9]), 'gx', label='Peak')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak2_log.pdf')
plt.clf()


k_B = 8.6173 * 10**(-5)


def gerade(t, a, b):
    return a * t + b


time = np.linspace(0, len(Temp_peak1), len(Temp_peak1))
params_temp_peak1 = ucurve_fit(gerade, time, Temp_peak1)
t_plot = np.linspace(0, 22, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_temp_peak1)), 'b-', label='Fit')
plt.plot(time, Temp_peak1, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_Time_Peak1.pdf')
plt.clf()

write('build/2_Heizrate_peak1.tex', make_SI(params_temp_peak1[0], r'\kelvin\per\minute', figures=1))


time = np.linspace(0, len(Temp_peak2), len(Temp_peak2))
params_temp_peak2 = ucurve_fit(gerade, time, Temp_peak2)
t_plot = np.linspace(0, 31, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_temp_peak2)), 'b-', label='Fit')
plt.plot(time, Temp_peak2, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_Time_Peak2.pdf')
plt.clf()


write('build/2_Heizrate_peak2.tex', make_SI(params_temp_peak2[0], r'\kelvin\per\minute', figures=1))


def gerade(invT, a, b):
    return (-a / k_B) * invT + b


params = ucurve_fit(gerade, 1 / (Temp_peak1[0:12]), noms(Strom_peak1_log[0:12]))   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(0.003854, 0.00415, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params)), 'b-', label='Fit')
plt.errorbar(1 / (Temp_peak1[0:12]), noms(Strom_peak1_log[0:12]), fmt='rx', yerr=stds(Strom_peak1_log[0:12]), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak_log_fit.pdf')
plt.clf()

tau_0 = (1 / params_temp_peak1[0]) * unp.exp(-params[0] / (k_B * Temp_peak1[12])) * k_B * Temp_peak1[12]**2 / params[0]
write('build/2_tau_peak1_variante1.tex', make_SI(tau_0 * 10**15, r'\pico\second', figures=1))       # type in Anz. signifikanter Stellen
write('build/2_tau_peak1_variante1_temp.tex', make_SI(Temp_peak1[12], r'\kelvin', figures=1))       # type in Anz. signifikanter Stellen

write('build/2_Energie_peak1_variante1.tex', make_SI(params[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen


params = ucurve_fit(gerade, 1 / (Temp_peak2[0:22]), noms(Strom_peak2_log[0:22]))   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(0.003175, 0.00355, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params)), 'b-', label='Fit')
plt.errorbar(1 / (Temp_peak2[0:22]), noms(Strom_peak2_log[0:22]), fmt='rx', yerr=stds(Strom_peak2_log[0:22]), label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [1/K]$')
plt.ylabel(r'$ln(I \:/\: [A])$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak2_log_fit.pdf')
plt.clf()


tau_0 = (1 / params_temp_peak2[0]) * unp.exp(-params[0] / (k_B * Temp_peak2[22])) * k_B * Temp_peak2[22]**2 / params[0]
write('build/2_tau_peak2_variante1.tex', make_SI(tau_0 * 10**15, r'\pico\second', figures=1))       # type in Anz. signifikanter Stellen
write('build/2_tau_peak2_variante1_temp.tex', make_SI(Temp_peak2[22], r'\kelvin', figures=1))       # type in Anz. signifikanter Stellen

write('build/2_Energie_peak2_variante1.tex', make_SI(params[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen

# zweites Verfahren


#Strom_peak1_integral = np.ones(len(Temp_peak1))
Strom_peak1_integral = unp.uarray(np.zeros(len(Temp_peak1)), np.zeros(len(Temp_peak1)))
Strom_peak2_integral = unp.uarray(np.zeros(len(Temp_peak2)), np.zeros(len(Temp_peak2)))

for j in range(1, len(Temp_peak1)):
    help = ufloat(0, 0)
    for i in range(j, len(Temp_peak1)):
        help = help + (Temp_peak1[i] - Temp_peak1[i - 1]) * (Strom_peak1[i] + Strom_peak1[i - 1]) / 2   # Trapezverfahren
    Strom_peak1_integral[j - 1] = help

for j in range(1, len(Temp_peak2)):
    help = ufloat(0, 0)
    for i in range(j, len(Temp_peak2)):
        help = help + (Temp_peak2[i] - Temp_peak2[i - 1]) * (Strom_peak2[i] + Strom_peak2[i - 1]) / 2   # Trapezverfahren
    Strom_peak2_integral[j - 1] = help


ln_Strom_peak1_mod = unp.log(Strom_peak1_integral[0:len(Strom_peak1_integral) - 1]) - unp.log(Strom_peak1[0:len(Strom_peak1) - 1])
ln_Strom_peak2_mod = unp.log(Strom_peak2_integral[0:len(Strom_peak2_integral) - 1]) - unp.log(Strom_peak2[0:len(Strom_peak2) - 1])


#plt.plot(Temp_peak1[0:15], ln_Strom_peak1_mod[0:15], 'bx', label='Peak')


grenze = len(ln_Strom_peak1_mod) - 1
plt.errorbar(Temp_peak1[0:len(Temp_peak1) - 1], noms(ln_Strom_peak1_mod), fmt='rx', yerr=stds(ln_Strom_peak1_mod), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$ln(I) \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak_log_2.pdf')
plt.clf()

plt.errorbar(Temp_peak2[0:len(Temp_peak2) - 1], noms(ln_Strom_peak2_mod), fmt='rx', yerr=stds(ln_Strom_peak2_mod), label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$ln(I) \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak2_log_2.pdf')
plt.clf()


def gerade(T, a, b):
    return (a / (T * k_B)) + b


params_peak1 = ucurve_fit(gerade, Temp_peak1[0:len(Temp_peak1) - 1], ln_Strom_peak1_mod)   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(241, 273, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_peak1)), 'b-', label='Fit')
plt.errorbar(Temp_peak1[0:len(Temp_peak1) - 1], noms(ln_Strom_peak1_mod), fmt='rx', yerr=stds(ln_Strom_peak1_mod), label='Peak1')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$ln(I) \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak_log_3.pdf')
plt.clf()

params_peak2 = ucurve_fit(gerade, Temp_peak2[0:len(Temp_peak2) - 1], ln_Strom_peak2_mod)   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(282, 328, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_peak2)), 'b-', label='Fit')
plt.errorbar(Temp_peak2[0:len(Temp_peak2) - 1], noms(ln_Strom_peak2_mod), fmt='rx', yerr=stds(ln_Strom_peak2_mod), label='Peak2')
plt.xlabel(r'$Temp \:/\: \: [K]$')
plt.ylabel(r'$ln(I) \:/\: [A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_current_peak2_log_3.pdf')
plt.clf()


def gerade(t, a, b):
    return a * t + b


time = np.linspace(0, len(Temp_peak1), len(Temp_peak1))
params_temp_peak1 = ucurve_fit(gerade, time, Temp_peak1)
t_plot = np.linspace(0, 22, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_temp_peak1)), 'b-', label='Fit')
plt.plot(time, Temp_peak1, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_Time_Peak1.pdf')
plt.clf()


time = np.linspace(0, len(Temp_peak2), len(Temp_peak2))
params_temp_peak2 = ucurve_fit(gerade, time, Temp_peak2)
t_plot = np.linspace(0, 31, 10)
plt.plot(t_plot, gerade(t_plot, *noms(params_temp_peak2)), 'b-', label='Fit')
plt.plot(time, Temp_peak2, 'rx', label='Messdaten')
plt.xlabel(r'$t \:/\: \: [min]$')
plt.ylabel(r'$Temp \:/\: [K]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/2_Temp_Time_Peak2.pdf')
plt.clf()

t_0 = (unp.exp(params_peak1[1]) / params_temp_peak1[0])
write('build/2_tau_peak1_variante2.tex', make_SI(t_0 * 10**18, r'\atto\second', figures=1))       # type in Anz. signifikanter Stellen
write('build/2_Energie_peak1_variante2.tex', make_SI(params_peak1[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen


t_0 = (unp.exp(params_peak2[1]) / params_temp_peak2[0])
write('build/2_tau_peak2_variante2.tex', make_SI(t_0 * 10**18, r'\atto\second', figures=1))       # type in Anz. signifikanter Stellen
write('build/2_Energie_peak2_variante2.tex', make_SI(params_peak2[0], r'\electronvolt', figures=2))       # type in Anz. signifikanter Stellen

########## IMPORT ###################################################################################################################################################
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
