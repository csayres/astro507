import numpy
numpy.seterr("ignore")
import matplotlib.pyplot as plt
from scipy.integrate import quad

kb = 1.3806504e-16 # erg per K
h = 6.62606885e-27 # erg per s
c = 2.99792458e10 # cm per s
jy = 10**-23 # erg s^-1 cm^-2 Hz^-1
T_CMB = 2.725 # K
data = numpy.loadtxt("firas_monopole_spec_v1.txt")
nu = data[:,0] * c # convert to s^-1
I = data[:,1] * (10**6*jy) # convert to specific intensity
sigmaI = data[:,3] * (10**3*jy) # convert to specific intensity

def blackbody(nu, mu, T):
    # nu is s^-1, frequency
    # mu is erg, chemical potential
    # T is K, temperature
    # returns units of erg s^-1 cm^-2 sr^-1 per frequency
    # specific intensity

    return 2*h*nu**3/c**2/(numpy.exp((h*nu-mu)/(kb*T)) - 1)
    # return 2*h*nu**3/c**2/(numpy.exp(h*nu/(kb*T)-mu/(kb*T)) - 1)


plt.errorbar(nu, I, yerr=sigmaI, alpha=0.5, color="b", label="data")
modeled = blackbody(nu, 0, T_CMB)
# plt.plot(nu, I, color="b")
plt.plot(nu, modeled, ".r", alpha=0.5, label="model")
plt.xlabel(r"$\nu$ $(s^{-1})$")
plt.ylabel(r"$I_\nu$")
plt.legend()
plt.savefig("bb.png")
plt.close()

# plot residual
plt.figure()
plt.errorbar(nu, I-modeled, yerr=sigmaI, label="data - model")
plt.xlabel(r"$\nu$ $(s^{-1})$")
plt.ylabel(r"$I_\nu$")
plt.legend()
plt.savefig("resid.png")
plt.close()

# construct a chi2 surface varying T and mu
tPoints = 500
muPoints = 300
Ts = numpy.linspace(T_CMB-1e-4, T_CMB+1e-4, tPoints)
muSet = []
chi2Grid = numpy.zeros((tPoints, muPoints))

bestChi = 1e16
for i,T in enumerate(Ts):
    Mus = numpy.linspace(-1e-4, 1e-4, muPoints)*kb*T
    muSet.append(Mus)
    for j,mu in enumerate(Mus):
        chi2 = numpy.sum((I - blackbody(nu, mu, T))**2/sigmaI**2)
        # reduced chi2
        if chi2 < bestChi:
            # print("bestChi", chi2)
            bestChi = chi2
        chi2Grid[i,j] = chi2
muSet = numpy.array(muSet)

chi2Grid = chi2Grid - numpy.min(chi2Grid)
bestFit = numpy.argwhere(chi2Grid==0).flatten()
# print("best Fit", bestFit)
chi2levels = [1, 2.71, 6.63]

fig = plt.figure()
im = plt.imshow(chi2Grid, origin="lower")
ax = plt.gca()
minChi = numpy.min(chi2Grid)
ax.contour(chi2Grid, levels=chi2levels, colors='k')#, origin='image')#, extent=extent)
plt.xlabel(r"$\mu k_BT$")
plt.ylabel(r"T (K)")
ax.set_xticks([0, 299])
ax.set_yticks([0, 499])
ax.set_xticklabels([r"(0 - 1e-4)$\mu k_b T$", r"(0 + 1e-4)$\mu k_b T$"])
ax.set_yticklabels([r"$T_{cmb}$ - 1e-4 K", r"$T_{cmb}$ + 1e-4 K"])
# fig.colorbar(im)
plt.plot(bestFit[1], bestFit[0], "Xr", label="best fit chi^2")
plt.legend()
plt.savefig("cont.png")
plt.close()

plt.figure(figsize=(10,5))
# extract the mus where chi2 < 6.63
goodMus = [muSet[x][y]/(kb*Ts[x]) for x,y in numpy.argwhere(chi2Grid<6.63)]
plt.hist(goodMus)
plt.xlabel(r"$\frac{\mu}{k_B T}$ values in 99 percent confidence interval")
plt.savefig("conf.png")
plt.close()


# problem 3
T = 10**6 # temp of star K
pc = 325 # g per cm^3 central density of star
me = 9.10938e-28 # g electron mass
mp = 1.67e-24 # g mass of proton
nh = pc/(mp+0.1*mp*4) # number density of hydrogen
nhe = 0.1*nh # number density of helium
ne = nh + 2*nhe # electron number density
print("ne", ne)
s = 1/2 # electron spin states
lambd = h/(2*numpy.pi*me*kb*T)**(0.5)

print("classical fugacity z: %.2f"%(ne*lambd**3))
# z fucacity = exp(mu/(kbT))

def fermiDirac(w, nu, z):
    return w**(nu)/(numpy.exp(w)/z + 1)

density = []
zs = numpy.logspace(-1,10,100000)
for z in zs:
    fz = quad(fermiDirac, 0, numpy.inf, args=(0.5, z))
    # get density
    _ne = 2*(2*s+1)/(numpy.pi**(1/2)*lambd**3)*fz[0]
    density.append(_ne)
density = numpy.array(density)

# find fugacity at given ne
minInd = numpy.argmin(numpy.abs(density-ne))
bestFugacity = zs[minInd]
print("solved for fucacity", bestFugacity)

plt.plot(numpy.log(zs), density, 'b', label="fermi-dirac")
plt.axhline(ne, color="r", label=r"$n_e$ = %.2e cm$^{-3}$"%(ne))
plt.plot(numpy.log(bestFugacity), ne, 'ok', markersize=10, label="log fugacity = %.2f"%(numpy.log(bestFugacity)))
plt.xlabel("log z")
plt.ylabel("number density")
plt.legend()
plt.savefig("3a.png")

# next find pressure
fz = quad(fermiDirac, 0, numpy.inf, args=(1.5, bestFugacity))
pressure = 4*(2*s+1)/(3*numpy.pi**(1/2)*lambd**3)*fz[0]*kb*T
print("pressure %.2e g cm^-1 s^-2"%pressure)

# are they relativistic?
# when fermi momentum is comparable to me*c
pf = h*(3*ne/(8*numpy.pi))**(1/3)
print("fermi momentum to me*c: %.2f << 1, so non-relativistic"%(pf/(me*c)))

# classical pressure
pClassical = ne*kb*T
# compare ratio:
print("ratio of classical pressure to degenerate pressure: %.2f"%(pClassical/pressure))

# estimate mass/radius
# R = (8.44*M/pc)**(1/3)
mu_e = pc/(mp*ne)
F = quad(fermiDirac, 0, numpy.inf, args=(0.5, bestFugacity))[0]
A = 0.13 / (mu_e*F)**(2/3)
M = (A**3*pc/8.44)**(1/2)
R = (8.44*M/pc)**(1/3)
print("Mass = %.2f solar masses"%M)
print("Radius = %.2f solor radii"%R)


