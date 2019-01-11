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

from io import StringIO
################################################ Finish importing custom libraries #################################################
U, I = np.genfromtxt('messdaten/IV_kurve.txt', unpack=True)

plt.plot(U, I, 'rx', label='Messdaten')
plt.xlabel(r'$U \: /\: [V]$')
plt.ylabel(r'$I \:/\: [\mu A]$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/IV_Kurve.pdf')
plt.clf()


f = open('messdaten/19_01_07_Kaiser_Hoetting/Pedestal.txt', 'r')
l = []
l = [line.split() for line in f]
# Pedestal = []
# for i in range(2, len(l), 1):
#     #Pedestal.append(np.genfromtxt(StringIO(l[i]), delimiter=";"))
ADC = []
for i in range(1, len(l)):
    y = []
    x = l[i][0].split(";")
    for item in x:
        y.append(float(item))
    ADC.append(y)

# print(len(ADC))
# print((Pedestal[0]))


# Argumentwerte als 1D Arrays erzeugen
x_1d = np.linspace(0, 127, 128)
y_1d = np.linspace(0, 999, 1000)

# # Argumentwerte als 2D Arrays erzeugen
y_2d, x_2d = np.meshgrid(y_1d, x_1d)

# # Interessante Daten erzeugen
# z_2d = 1/(x_2d**2 + y_2d**2 + 1) * np.cos(np.pi * x_2d) * np.cos(np.pi * y_2d)

# plt.figure()
# plt.pcolormesh(x_2d, y_2d, ADC, cmap=plt.get_cmap("plasma"))
# plt.colorbar()
# plt.xlabel('Channels')
# plt.ylabel('Events')
# # plt.gca().set_aspect("equal")
# plt.savefig('build/ADC.pdf')
# plt.clf()

Pedestals = []
for i in range(len(ADC)):
    help = 0
    for j in range(len(ADC[i])):
        help = help + ADC[i][j]
    Pedestals.append((1 / len(ADC[i])) * help)

Common_Mode_Shift = []
for i in range(len(ADC[0])):
    help = 0
    for j in range(len(ADC)):
        help = help + (ADC[j][i] - Pedestals[j])
    Common_Mode_Shift.append((1 / 128) * help)

Noise = []
for i in range(len(ADC)):
    help = 0
    for j in range(len(ADC[i])):
        help = help + (ADC[i][j] - Pedestals[i] - Common_Mode_Shift[j])**2
    Noise.append(np.sqrt((1 / (len(ADC[i]) - 1)) * help))


# plt.bar(x_1d, Noise, align='center')
# plt.xlabel('Channels')
# plt.ylabel('Noise')
# # plt.gca().set_aspect("equal")
# plt.savefig('build/Noise.pdf')
# plt.clf()

# plt.bar(x_1d, Pedestals, align='center')
# plt.xlabel('Channels')
# plt.ylabel('Pedestals')
# # plt.gca().set_aspect("equal")
# plt.savefig('build/Pedestals.pdf')
# plt.clf()

# plt.bar(y_1d, Common_Mode_Shift, align='center')
# plt.xlabel('Events')
# plt.ylabel('CommonModeShift')
# # plt.gca().set_aspect("equal")
# plt.savefig('build/Common_Mode_Shift.pdf')
# plt.clf()

Cluster = np.linspace(-10, 11, 85)
Cluster_Result = np.zeros(len(Cluster))
for i in range(len(Common_Mode_Shift)):
    for j in range(len(Cluster) - 1):
        if Cluster[j] < Common_Mode_Shift[i] and Cluster[j + 1] > Common_Mode_Shift[i]:
            Cluster_Result[j] = Cluster_Result[j] + 1
sum = 0
for i in range(len(Cluster_Result)):
    sum = sum + Cluster_Result[i]
for i in range(len(Cluster_Result)):
    Cluster_Result[i] = Cluster_Result[i] / sum

plt.bar(Cluster, Cluster_Result, align='center')
plt.xlabel('Häugikeit')
plt.ylabel('CommonModeShift')
# plt.gca().set_aspect("equal")
plt.savefig('build/CMS_Gauß.pdf')
plt.clf()


puls, ADC_CH10 = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/Calib/Kalibration_b_Ch10_90.txt', skip_header=2, delimiter=' ', unpack=True)
puls, ADC_CH10_0Volt = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/Calib/Kalibration_b_Ch10.txt', skip_header=2, delimiter=' ', unpack=True)
puls, ADC_CH35 = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/Calib/Kalibration_b_Ch35.txt', skip_header=2, delimiter=' ', unpack=True)
puls, ADC_CH70 = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/Calib/Kalibration_b_Ch70.txt', skip_header=2, delimiter=' ', unpack=True)
puls, ADC_CH100 = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/Calib/Kalibration_b_Ch100.txt', skip_header=2, delimiter=' ', unpack=True)
puls, ADC_CH120 = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/Calib/Kalibration_b_Ch120.txt', skip_header=2, delimiter=' ', unpack=True)

ADC_mean = []
for i in range(210):
    sum = 0
    sum = (ADC_CH10[i] + ADC_CH35[i] + ADC_CH70[i] + ADC_CH100[i] + ADC_CH120[i]) / 5
    ADC_mean.append(sum)


plt.plot(puls[0:210], ADC_CH10[0:210], 'b.', label='Ch10')
plt.plot(puls[0:210], ADC_CH10_0Volt[0:210], 'r-', label='Ch10 0V')
plt.plot(puls[0:210], ADC_CH35[0:210], 'g.', label='Ch35')
plt.plot(puls[0:210], ADC_CH70[0:210], 'm.', label='Ch70')
plt.plot(puls[0:210], ADC_CH100[0:210], 'y.', label='Ch100')
plt.plot(puls[0:210], ADC_CH120[0:210], 'c.', label='Ch120')
plt.plot(puls[0:210], ADC_mean[0:210], 'kx', label='Mean')
plt.xlabel('Puls in e')
plt.ylabel('ADC')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Calib.pdf')
plt.clf()


plt.plot(puls[50:100], ADC_CH10[50:100], 'b.', label='Ch10')
plt.plot(puls[50:100], ADC_CH10_0Volt[50:100], 'r-', label='Ch10 0V')
plt.plot(puls[50:100], ADC_CH35[50:100], 'g.', label='Ch35')
plt.plot(puls[50:100], ADC_CH70[50:100], 'm.', label='Ch70')
plt.plot(puls[50:100], ADC_CH100[50:100], 'y.', label='Ch100')
plt.plot(puls[50:100], ADC_CH120[50:100], 'c.', label='Ch120')
plt.plot(puls[50:100], ADC_mean[50:100], 'kx', label='Mean')
plt.xlabel('Puls in e')
plt.ylabel('ADC')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Calib_detail.pdf')
plt.clf()


def f(x, a, b, c, d, e):
    return a * x**4 + b * x**3 + c * x**2 + d * x + e


# for i in range(len(puls)):
#     puls[i] = puls[i] / 1000

plt.plot(ADC_mean[0:210], puls[0:210], 'kx', label='Mean')
params = ucurve_fit(f, ADC_mean[0:150], puls[0:150], p0=[0.001, -0.03, 1, 460, 0])
# p0=[0.001, -0.03, 1, 460, 0]   # p0 bezeichnet die Startwerte der zu fittenden Parameter
t_plot = np.linspace(0, 300, 100)
plt.plot(t_plot, f(t_plot, *noms(params)), 'b-', label='Fit')
plt.axvline(x=250, ymin=0, ymax=200, color='r', linestyle='--', label='Grenze')
#plt.axvline(x=150000, ymin=0, ymax=270, color='r', linestyle='--', label='Grenze')
plt.xlabel('ADC')
plt.ylabel('Puls in ke')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Calib_fit.pdf')
plt.clf()
Params_Umrechnung = params
# a = (0.000 030 ± 0.000 002) e
# b = (−0.0084 ± 0.0008) e
# c = (1.1 ± 0.1) e
# d = (262 ± 8) e
# e = (25 ± 137) e


f = open('messdaten/19_01_07_Kaiser_Hoetting/Laserscan.txt', 'r')
l = []
l = [line.split() for line in f]
ADC = []
for i in range(1, len(l)):
    y = []
    x = l[i][0].split(",")
    for item in x:
        y.append(float(item))
    ADC.append(y)

# print(len(ADC))
# print(len(ADC[1]))


# Argumentwerte als 1D Arrays erzeugen
x_1d = np.linspace(0, 34, 35)
y_1d = np.linspace(0, 127, 128)

# # Argumentwerte als 2D Arrays erzeugen
y_2d, x_2d = np.meshgrid(y_1d, x_1d)

# # Interessante Daten erzeugen
# z_2d = 1/(x_2d**2 + y_2d**2 + 1) * np.cos(np.pi * x_2d) * np.cos(np.pi * y_2d)

plt.figure()
plt.pcolormesh(x_2d, y_2d, ADC, cmap=plt.get_cmap("plasma"))
plt.colorbar()
plt.xlabel('Abstand')
plt.ylabel('Channels')
# plt.gca().set_aspect("equal")
plt.savefig('build/Laserscan_komplett.pdf')
plt.clf()

ADC_zoom4555 = []
for i in range(len(ADC)):
    help = []
    for j in range(45, 56):
        help.append(ADC[i][j])
    ADC_zoom4555.append(help)

x_1d = np.linspace(0, 34, 35)
y_1d = np.linspace(45, 55, 11)

# # Argumentwerte als 2D Arrays erzeugen
y_2d, x_2d = np.meshgrid(y_1d, x_1d)

# # Interessante Daten erzeugen
# z_2d = 1/(x_2d**2 + y_2d**2 + 1) * np.cos(np.pi * x_2d) * np.cos(np.pi * y_2d)

plt.figure()
plt.pcolormesh(x_2d, y_2d, ADC_zoom4555, cmap=plt.get_cmap("plasma"))
plt.colorbar()
plt.xlabel('Abstand')
plt.ylabel('Channels')
# plt.gca().set_aspect("equal")
plt.savefig('build/Laserscan_zoom.pdf')
plt.clf()


ADC_zoom4953 = []
for i in range(49, 52):
    help = []
    for j in range(len(ADC)):
        help.append(ADC[j][i])
    ADC_zoom4953.append(help)
# print(len(ADC_zoom4953))
# print(len(ADC_zoom4953[0]))

pos = np.linspace(0, 34, 35)
plt.plot(pos, ADC_zoom4953[0], 'r-', label='Channel 50')
plt.plot(pos, ADC_zoom4953[1], 'b-', label='Channel 51')
plt.plot(pos, ADC_zoom4953[2], 'g-', label='Channel 52')
plt.axvline(x=4, ymin=0, ymax=160, color='c', linestyle='--', label='Streifenbreite')
plt.axvline(x=10.95, ymin=0, ymax=160, color='c', linestyle='--')
plt.axvline(x=11.05, ymin=0, ymax=160, color='m', linestyle='--', label='Streifenabstand')
plt.axvline(x=19.95, ymin=0, ymax=160, color='m', linestyle='--')
plt.axvline(x=20.05, ymin=0, ymax=160, color='c', linestyle='--', label='Streifenbreite')
plt.axvline(x=27, ymin=0, ymax=160, color='c', linestyle='--')
plt.xlabel(r'Position$\:/\: [10^{-5} m]$')
plt.ylabel('ADC')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Laserscan_Pos.pdf')
plt.clf()

##########################################################################################################################################################################

CCE_read = np.linspace(0, 200, 21)
CCE = []
for i in CCE_read:
    CCE.append(np.genfromtxt("messdaten/19_01_07_Kaiser_Hoetting/CCEL/" + str(int(i)) + "CCEL.txt", skip_header=1, unpack=True))

Channel = np.linspace(48, 53, 6)
c = []  # color of points
s = []

for i in range(len(CCE_read)):
    s.append(305 - CCE_read[i] * 1.5)
    c.append(CCE_read[i] * 0.05)

c = []
x_48 = []
x_49 = []
x_50 = []
x_51 = []
x_52 = []
signal_48 = []
signal_49 = []
signal_50 = []
signal_51 = []
signal_52 = []
# print('HIER')
# print((len(CCE)))
for i in range(len(CCE)):
    signal_48.append(CCE[i][48])
    signal_49.append(CCE[i][49])
    signal_50.append(CCE[i][50])
    signal_51.append(CCE[i][51])
    signal_52.append(CCE[i][52])
    x_48.append(int(48))
    x_49.append(int(49))
    x_50.append(int(50))
    x_51.append(int(51))
    x_52.append(int(52))
    c.append(i)


cm = plt.cm.get_cmap('rainbow')
plt.scatter(x_50, signal_50, c=c, s=s, cmap=plt.cm.rainbow)
plt.scatter(x_51, signal_51, c=c, s=s, cmap=plt.cm.rainbow)
plt.scatter(x_52, signal_52, c=c, s=s, cmap=plt.cm.rainbow)
plt.scatter(x_48, signal_48, c=c, s=s, cmap=plt.cm.rainbow)
sc = plt.scatter(x_49, signal_49, c=c, s=s, cmap=cm)
clb = plt.colorbar(sc)
clb.set_label('label', labelpad=13, y=0.5, rotation=270)
plt.savefig('build/test.pdf')
plt.clf()


maximum = max(signal_50)
for i in range(len(signal_50)):
    signal_50[i] = signal_50[i] / maximum


# def CCE_function(U, U_dep, a):
#     D = 300 * 10**(-6)
#     # U_dep = 60
#     return ((1 - np.exp(-D * np.sqrt((U / U_dep)) / a)) / (1 - np.exp(-D / a)))

def CCE_function(U, U_depl, a):  # a: mitt. Eindringtiefe
    dicke = 300 * 1e-6  # Sensordicke 300 micrometer
    depl_zone = np.array([dicke * min([np.sqrt(u / U_depl), 1]) for u in U])  # Dicke der Depletionszone (Gleichung 1 im Protokoll)
    func = (1 - np.exp(-depl_zone / a)) / (1 - np.exp(-dicke / a))  # Gleichung 2 im Protokoll
    return func


b = ((0, 0), (200, 300 * 1e-6))
params = ucurve_fit(CCE_function, CCE_read, signal_50, bounds=b)   # p0 bezeichnet die Startwerte der zu fittenden Parameter
U_plot = np.linspace(0, 200, 100)
plt.plot(U_plot, CCE_function(U_plot, *noms(params)), 'b-', label='Fit')

# print(params)

plt.plot(CCE_read, signal_50, 'kx', label='Mean')
plt.xlabel(r'$U \: /\: [V]$')
plt.ylabel('Intensität')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/CCE_fit.pdf')
plt.clf()


# f = open('messdaten/19_01_07_Kaiser_Hoetting/CCEQ/10_Cluster_adc_entries.txt', 'r')
# l = []
# l = [line.split() for line in f]
# ADC = []
# for i in range(1, len(l)):
#     y = []
#     for item in l[i]:
#         y.append(float(item))
#     ADC.append(y)


def Funktion_Auslese(Datei_Pfad):
    f = open(Datei_Pfad, 'r')
    l = []
    l = [line.split() for line in f]
    help = []
    for i in range(1, len(l)):
        for item in l[i]:
            help.append(float(item))
    return help


volt = np.linspace(0, 200, 21)
pfad_begin = 'messdaten/19_01_07_Kaiser_Hoetting/CCEQ/'
pfad_end = '_Cluster_adc_entries.txt'
CCEQ = []
for i in volt:
    CCEQ.append(Funktion_Auslese(pfad_begin + str(int(i)) + pfad_end))
# print((CCEQ[0][85:95]))

CCEQ_mean = []
for i in range(len(CCEQ)):
    CCEQ_mean.append(np.mean(CCEQ[i]))

maximum = max(CCEQ_mean)
for i in range(len(CCEQ_mean)):
    CCEQ_mean[i] = CCEQ_mean[i] / maximum


b = ((0, 0), (200, 300 * 1e-6))
params = ucurve_fit(CCE_function, volt, CCEQ_mean, bounds=b)   # p0 bezeichnet die Startwerte der zu fittenden Parameter
U_plot = np.linspace(0, 200, 100)
plt.plot(U_plot, CCE_function(U_plot, *noms(params)), 'b-', label='Fit')

# print(params)

plt.plot(volt, CCEQ_mean, 'kx', label='Mean')
plt.xlabel(r'$U \: /\: [V]$')
plt.ylabel('Intensität')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/CCEQ_fit.pdf')
plt.clf()


frequency = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/cluster_size.txt', skip_header=1, unpack=True)
channel = np.linspace(0, 127, 128)
total = 0
for i in range(len(frequency)):
    total = total + frequency[i]

for i in range(len(frequency)):
    frequency[i] = frequency[i] / total
# plt.bar(x_1d, Noise, align='center')
plt.bar(channel[0:16], frequency[0:16])
plt.xlabel('Anzahl der Kanäle pro Event')
plt.ylabel('normierte Häufigkeit')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Quellenmessung_frequeny_channel.pdf')
plt.clf()


frequency = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/number_of_clusters.txt', skip_header=1, unpack=True)
channel = np.linspace(0, 127, 128)
total = 0
for i in range(len(frequency)):
    total = total + frequency[i]

for i in range(len(frequency)):
    frequency[i] = frequency[i] / total
# plt.bar(x_1d, Noise, align='center')
plt.bar(channel[0:6], frequency[0:6])
plt.xlabel('Anzahl der Cluster pro Event')
plt.ylabel('normierte Häufigkeit')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Quellenmessung_frequeny_cluster.pdf')
plt.clf()


frequency = np.genfromtxt('messdaten/19_01_07_Kaiser_Hoetting/hitmap.txt', skip_header=1, unpack=True)
channel = np.linspace(0, 127, 128)
total = 0
for i in range(len(frequency)):
    total = total + frequency[i]

for i in range(len(frequency)):
    frequency[i] = frequency[i] / total
# plt.bar(x_1d, Noise, align='center')
plt.bar(channel, frequency)
plt.xlabel('Kanäle')
plt.ylabel('normierte Anzahl der Hits')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Quellenmessung_hitmap.pdf')
plt.clf()


def function(x, a, b, c, d, e):
    return a * x**4 + b * x**3 + c * x**2 + d * x + e


f = open('messdaten/19_01_07_Kaiser_Hoetting/Cluster_adc_entries.txt', 'r')
l = []
l = [line.split() for line in f]
ADC = []
for i in range(1, len(l)):
    y = []
    for item in l[i]:
        y.append(float(item))
    ADC.append(y)


for i in range(len(ADC)):
    x = 0
    for j in range(len(ADC[i])):
        ADC[i][j] = 3.6 * (function(ADC[i][j], *noms(Params_Umrechnung)))


ADC_sum = []
for i in range(len(ADC)):
    x = 0
    for j in range(len(ADC[i])):
        x = x + ADC[i][j]
    ADC_sum.append(x)

Energie_mean = np.mean(ADC_sum)
Energie_max = max(ADC_sum)

print(Energie_max)
print(Energie_mean)
Cluster = np.linspace(0, 1215230, 400000)
Cluster_Result = np.zeros(len(Cluster))
ADC_sum.sort()
print(ADC_sum[0:10])
# for i in range(len(ADC_sum)):
#     for j in range(len(Cluster) - 1):
#         if Cluster[j] < ADC_sum[i] and Cluster[j + 1] > ADC_sum[i]:
#             Cluster_Result[j] = Cluster_Result[j] + 1
#             break

# for i in range(len(Cluster)):
#     Cluster[i] = Cluster[i] / 1000
# item = 0
# Cluster_Result_max = max(Cluster_Result)
# for i in range(len(Cluster_Result)):
#     if Cluster_Result[i] == Cluster_Result_max:
#         item = i
#         break

# plt.axvline(x=Energie_mean / 1000, ymin=0, ymax=9000, color='c', linestyle='--', label='Mittelwert')
# plt.axvline(x=Cluster[i], ymin=0, ymax=9000, color='m', linestyle='--', label='MPV')
# plt.axvline(x=function(250,*noms(Params_Umrechnung)), ymin=0, ymax=9000, color='r', linestyle='--', label='Grenzwert')
# plt.plot(Cluster, Cluster_Result, 'bx')
# plt.xlabel('Energie pro Cluster in keV')
# plt.ylabel('Häufigkeit')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/Quellenmessung_Cluster.pdf')
# plt.clf()
# ################################ FREQUENTLY USED CODE ################################
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
