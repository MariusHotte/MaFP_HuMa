import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as op
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
plt.rcParams['figure.figsize']= (10,8)
plt.rcParams['font.size']=16
plt.rcParams['lines.linewidth']=0.5

x1, y1= np.genfromtxt('Messwerte/Stromkurve1.txt', unpack= True)
x1= x1*100
x2, y2= np.genfromtxt('Messwerte/Stromkurve2.txt', unpack= True)
x2= x2*100

def func(x,a,b):
	return a*x+b
def func2(x,a,c):
	return a*x**2+c
def plot0(x1,y1,x2,y2, k, xlab, ylab):
	plt.figure()
	plt.plot(x1,y1,'o')
	plt.plot(x2,y2,'r')
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.grid()
	plt.savefig('Bilder/plot'+k+'.png')
def plot1(x1,y1,y1f,x2,y2, k,xlab,ylab):
	plt.figure()
	plt.errorbar(x1, y1,yerr= y1f, fmt='o')
	plt.plot(x2,y2,'r')
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.ylim(0,0.2e-22)
	plt.grid()
	plt.savefig('Bilder/plot'+k+'.png')
def plot2(x1,y1,y1f, k):
	plt.figure()
	plt.errorbar(x1, y1,yerr= y1f, fmt='o')
	plt.xlabel(r"$\Theta$")
	plt.ylabel(r"$N$")
	plt.grid()
	plt.savefig('Bilder/plot'+k+'.png')
def plot3(x1,y1,k):
	plt.figure()
	plt.plot(x1,y1,'o')
	plt.xlabel(r"")
	plt.ylabel(r"")
	plt.grid()
	plt.savefig('Bilder/plot'+k+'.png')
def plot4(x1,y1,y1f,x2,y2, k,xlab,ylab):
	plt.figure()
	plt.errorbar(x1, y1,yerr= y1f, fmt='o')
	plt.plot(x2,y2,'r')
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.grid()
	plt.savefig('Bilder/plot'+k+'.png')

params1, covariance1 = curve_fit(func, x1, y1)
errors1 = np.sqrt(np.diag(covariance1))
t= np.linspace(x1.min(),x1.max())
plot0(x1,y1,t,func(t,*params1),'1',r"$p$ in bar",r"$U$ in V")
params2, covariance2 = curve_fit(func, x2, y2)
errors2 = np.sqrt(np.diag(covariance2))
plot0(x2,y2,t,func(t,*params2),'2',r"$p$ in bar",r"$U$ in V")

x3, y3= np.genfromtxt('Messwerte/Diff.txt', unpack= True)
y3= y3
z3= np.sqrt(y3)

eps_0= 8.85e-12
e0= 1.602e-19
Z= 79
z= 2
E= 5.638e6*e0
Theta= np.linspace(x3.min(),x3.max(),1000)
ds= ((z*Z*e0**2)/(16*np.pi*eps_0*E))**2/(np.sin(Theta/2/180*np.pi))**4
Nd_0= 2022
n= 5.8e22*2.5
dN= y3*4/(n*Nd_0)
dNf= z3*4/(n*Nd_0)

plot2(x3,y3,z3,'3')
plot1(x3,dN,dNf,Theta, ds,'4',r"$\Theta$",r"$\frac{d\sigma}{d\Omega}$")

n2= n/2.5*4
d4N0= 963*4/(n2*Nd_0)
d4N7= 766*4/(n2*Nd_0)


X= [13,79,83]
nb= (9.78*1000)/(208.9*1.6e-27)
nal= (2.7*1000)/(26.9*1.6e-27)
Y= [546/400/(nal*3e-6),614/100/n, 603/400/(nb*1e-6)]
Yf= [np.sqrt(546)/400/(nal*3e-6),np.sqrt(614)/100/n, np.sqrt(603)/400/(nb*1e-6)]

params3, covariance3 = curve_fit(func2, X, Y)
errors3 = np.sqrt(np.diag(covariance3))
print(errors3)
t2= np.linspace(10,90)
plot4(X,Y,Yf,t2,func2(t2,*params3),'5',r"$Z$",r"$I_\alpha/(N \cdot a)$")