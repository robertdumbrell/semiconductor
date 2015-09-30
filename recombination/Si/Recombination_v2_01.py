# -*- coding: utf-8 -*-

# what this file should do:
#   define default parameters (which can be changed) for a wafer and recommbination processes
#   SRH, Auger, Radiative, Jo_emitter
#   plot tau(deltan) from model
#   fit data with any/all tau models

# CHANGE LOG:
# v2.01
#   The plotall functions now work
#   made dummy dependance of Deltan, so fitting works but functions call self.Deltan
#       No dummy Delta was left for radiative or Auger as no fitting of them is required
#
#   Constants Class
#           Removed def int, and self. before all items in constants class
#               This lets you be able to change constants and it affect things downstream, else it is not able to change things like doping - I would like to discuss this change, Mattias
#
#   Fitting Class
#           Fitting Graphs corrected
#           vairables Y1 and Y2 changed to tau_SHR and itau_SHR
#           Made function only accept positive lifetime values, is negtive lifetime possible?
#           Given inital values for Fitting, to help with convergance
#   Radiative Recombination has been updated to work properly.


# v2
#   inheritance:
#       Grandaddy = Constants()
#       Parents = SRH(), Auger(), Radiative(), Emitter()
#       Child = Lifetime()
#
#   changed "Constants()." references to "self."
#
#   "ReturnEqulibriumAndCarrierDensity".... really?
#   "vairables"...wut
#   changed SHR to SRH
#
#   don't need to include class "attributes" in "method" argument list
#
#   names of "Lifetime" methods changed for ambiguity reasons


# this is not polished yet, but at least it works, kinda


import sys
import re
import numpy
from numpy import *
from matplotlib.pylab import *

#from scipy.optimize import curve_fit
import scipy.optimize


class Constants(object):

    # def __init__(self):
    # constantly constant contants, lol
    kb = 1.3806488e-23  # J/K
    q = 1.6e-19  # C
    # not so constant constants
    T = 300  # K
    ni = 9.65e9  # /cm3 from memory need reference, make a function of T?
    Doping = 5e15  # /cm3
    Deltan = logspace(11, 19, 1000)
    MajorityCarrier = 'p'
    Width = 0.018
    Eg = 1.12  # eV
    # cm/s Thermal velocity from W. M. Bullis and H. R. Huff, J. Electrochem.
    # Soc. 143, 1399 ͑1996͒.
    vth = 1.1e7

    def Vth(self):
        return self.kb * self.T / self.q

    def n_and_p(self, Deltan):
        if (self.MajorityCarrier == 'n'):
            n = self.Doping + Deltan
            p = Deltan
            n_0 = self.Doping
            p_0 = self.ni**2 / self.Doping
        elif(self.MajorityCarrier == 'p'):
            p = self.Doping + Deltan
            n = Deltan
            p_0 = self.Doping
            n_0 = self.ni**2 / self.Doping
        return p, n, p_0, n_0


class SRH(Constants):

    def __init__(self):
        # default parameters
        Constants.__init__(self)
        self.Et = 0  # eV  defined as Et = E_trap - E_i
        self.tau_n = 1e-5  # s
        self.tau_p = 1e-4  # s

    def Cal_tau(self, N):
        # print N,self.sigma_n,self.vth, self.Et
        self.tau_n = 1. / N / self.sigma_n / self.vth  # s
        self.tau_p = 1. / N / self.sigma_p / self.vth  # s

    def Iron(self, N):
        # Taken from Iron detection in crystalline silicon by carrier lifetime
        # measurements for arbitrary injection and doping 2004
        self.Et = -self.Eg / 2 + .38  # eV  defined as Et = E_trap - E_i
        self.sigma_n = 5e-14  # s
        self.sigma_p = 7e-17  # s
        self.Cal_tau(N)

    def IronBoron(self, N):
        # Taken from Iron detection in crystalline silicon by carrier lifetime
        # measurements for arbitrary injection and doping 2004
        self.Et = self.Eg / 2 - .23  # eV  defined as Et = E_trap - E_i
        self.sigma_n = 3e-14  # s
        self.sigma_p = 2e-15  # s
        self.Cal_tau(N)

    def tau_SRH(self, Deltan=[0, 0], tau_n=0, tau_p=0, Et=0):
        if all(Deltan != 0):
            self.Deltan = Deltan
            # print 'inputed Deltan values used'
        if(tau_n != 0):
            self.tau_n = tau_n
            # print 'inputed tau_n values used',
        if(tau_p != 0):
            self.tau_p = tau_p
            # print 'inputed tau_p values used'
        # print self.tau_n,self.tau_p,self.Et

        return self.SchroderTextbook_SRH(self.tau_n, self.tau_p, self.Et)

    def SchroderTextbook_SRH(self, tau_n, tau_p, Et=0):
        # this is been confirmed to agree with SZE as well as PVlighthouse
        # It is assumed that Nt << p, n, and that the semiconductor is not
        # degenerate.
        p1 = self.ni * numpy.exp(-(Et) / self.Vth())
        n1 = self.ni * numpy.exp((Et) / self.Vth())

        p, n, p_0, n_0 = self.n_and_p(self.Deltan)

        return (tau_p * (n + n1) + tau_n * (p + p1)) / (p)

    def itau_SRH(self, Deltan=0, tau_n=0, tau_p=0, Et=0):
        # print tau_n,tau_p,Et
        return 1. / self.tau_SRH(Deltan, tau_n, tau_p, Et)

    def PlotAll(self, plt):

        ax = plt.add_subplot(111)
        ax.plot(inf, inf, 'k-.', label='SRH')
        ax.plot(self.Deltan, self.SchroderTextbook_SRH(
            self.tau_n, self.tau_p, self.Et), '-.', label='Schroder Textbook')
        legend(loc=0)
        loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('Lifetime (s)')
        # show()


class Emitter_Recombination(Constants):
    # To use this you need to input a J0, taken from Cuevas2004
    # this is also for a symetrical sample

    J0e = 1e-13

    def __init__(self):
        Constants.__init__(self)

    def Cuevas2004(self):
        return 2 * self.J0e / self.q / self.ni**2 / self.Width * (self.Deltan + self.Doping)

    def tau_EmitterRecombination(self):
        return 1. / self.Cuevas2004()

    def itau_EmitterRecombination(self):
        return 1. / self.tau_EmitterRecombination()

    def PlotAll(self, plt):

        ax = plt.add_subplot(111)

        ax.plot(self.Deltan, self.tau_EmitterRecombination(), ':',
                label='Emitter Recombination:' + str(self.J0e))

        legend(loc=0)
        loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('$\Tau_eff$ (s)')
        # show()


class Auger(Constants):

    def __init__(self):
        Constants.__init__(self)

    def tau_aug(self):
        return self.Richter2012()

    def itau_aug(self):
        return 1 / self.tau_aug()

    def SchroderTextbook(self):
        C1 = 1.8e-24
        C2 = 6e-25
        C3 = 3e-27
        return 1 / ((self.Doping + self.Deltan) *
                    (C2 * self.Doping**(0.65) + C3 * self.Deltan**(0.8)))

    def Richter2012(self):

        K_n,        K_p,    K_Delta, delta,  L_eeh,  N_eeh,  K_eeh,  P_ehh,  L_ehh,  K_ehh\
            = 2.5e-31,   8.5e-32, 3e-29,  .92,    13,     .66,    3.3e17, .63,    7.5,    7e17
        p, n, p_0, n_0 = self.n_and_p(self.Deltan)

        g_eeh = 1 + L_eeh * (1 - tanh((n_0 / K_eeh)**(N_eeh)))
        g_ehh = 1 + L_ehh * (1 - tanh((p_0 / K_ehh)**(P_ehh)))

        return self.Deltan / (n * p - self.ni**2) /\
            (K_n * g_eeh * n_0 + K_p * g_ehh *
             p_0 + K_Delta * self.Deltan**delta)

    def KerrCuevas2002(self):
        '''Incorrect delta value '''
        K_n, K_p, K_Delta, delta, L_eeh, N_eeh, K_eeh, P_ehh, L_ehh, K_ehh\
            = 2.8e-31, 9.9e-32, 3.79e-29, .92, 44, .29, 2e16, .29, 44, 2e16
        p, n, p_0, n_0 = self.n_and_p(self.Deltan)
        g_eeh = 1 + L_eeh * (1 - tanh((n_0 / K_eeh)**(N_eeh)))
        g_ehh = 1 + L_ehh * (1 - tanh((p_0 / K_ehh)**(P_ehh)))

        return self.Deltan / (n * p - self.ni**2) /\
            (K_n * g_eeh * n_0 + K_p * g_ehh *
             p_0 + K_Delta * self.Deltan**delta)

    def PlotAll(self, plt):

        ax = plt.add_subplot(111)
        # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
        ax.plot(self.Deltan, self.SchroderTextbook(), label='SchroderTextbook')
        ax.plot(self.Deltan, self.KerrCuevas2002(), label='KerrCuevas 2002')
        ax.plot(self.Deltan, self.Richter2012(), label='Richter 2012')

        ax.legend(loc=0)
        loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('Lifetime (s)')
        # show()


class Radiative(Constants):
    Blow = 4.73E-15

    def __init__(self):
        Constants.__init__(self)

    def tau_rad(self):
        return self.Altermatt2005()

    def B(self, deltan, Nd, T=300):

        self.T = T
        bmax, rmax, smax, wmax = 1.,     0.2,  1.5e18, 4e18
        rmin, smin, wmin = 0.,     1e7,   1e9
        b2, r1, s1, w1 = 0.54,   320.,  550.,  365.
        b4, r2, s2, w2 = 1.25,    2.5,  3.,     3.54

        bmin = rmax + (rmin - rmax) / (1. + (self.T / r1)**r2)
        b1 = (smax + (smin - smax) / (1. + (self.T / s1)**s2)) * 2
        b3 = (wmax + (wmin - wmax) / (1. + (self.T / w1)**w2)) * 2

        # print bmin

        return self.Blow * (bmin + (bmax - bmin) / (1 + ((2 * deltan + Nd) / b1)**b2 + ((2 * deltan + Nd) / b3)**b4))

    def Altermatt2005(self):
        # The version appearing in xxx, belive this one is correct. Though b1
        # and B3 may need to divded by 2.

        # bmax,rmax,smax,wmax=        1.,     0.2,  1.5e18, 4e18
        # rmin,smin,wmin=                     0.,     1e7,   1e9
        # b2,r1,s1,w1 =               0.54,   320.,  550.,  365.
        # b4,r2,s2,w2 =               1.25,    2.5,  3.,     3.54

        # bmin = rmax + (rmin-rmax)/(1.+(self.T/r1)**r2)
        # b1 = (smax + (smin-smax)/(1.+(self.T/s1)**s2))*2
        # b3 = (wmax + (wmin-wmax)/(1.+(self.T/w1)**w2))*2
        # Blow = 4.73E-15

        # This is a radiative coefficients

        # p,n,p_0,n_0 = self.n_and_p(self.Deltan)

        # B = (bmin+(bmax-bmin)/( 1+ ((2*self.Deltan+self.Doping)/b1)**b2+\
        # ((2*self.Deltan+self.Doping)/b3)**b4))\
        # *Blow

        # B = (bmin+(bmax-bmin)/( 1+ ((n+p)/b1)**b2+\
        #                             ((n+p)/b3)**b4))\
        #     *Blow
        # print B[0]
        return self.Deltan / (n * p - self.ni**2) / B(self.Deltan, self.Doping, self.T)

    def PlotAll(self, plt):

        ax = plt.add_subplot(111)

        ax.plot(inf, inf, 'k--', label='Radiative')
        ax.plot(self.Deltan, self.Altermatt2005(), '--', label='Altermatt2005')

        legend(loc=0)
        loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('Lifetime (s)')
        # show()

    def itau_rad(self):
        return 1 / self.tau_rad()


class TotalLifetime(Auger, SRH, Emitter_Recombination, Radiative):

    def __init__(self):
        # default parameters
        Constants.__init__(self)
        self.Et = 0  # eV  defined as Et = E_trap - E_i
        self.tau_n = 1e-5  # s
        self.tau_p = 1e-4  # s

    def PlotiComponants(self, plt=False):
        if plt == False:
            plot(self.Deltan, self.itau_SRH(), '--', label='SRH')
            plot(self.Deltan, self.itau_aug(), '--', label='Auger')
            plot(self.Deltan, self.itau_rad(), '--', label='Radiative')
            plot(self.Deltan, self.itau_EmitterRecombination(),
                 '--', label='Emitter')
        else:
            plt.plot(self.Deltan, self.itau_SRH(), '--', label='SRH')
            plt.plot(self.Deltan, self.itau_aug(), '--', label='Auger')
            plt.plot(self.Deltan, self.itau_rad(), '--', label='Radiative')
            plt.plot(
                self.Deltan, self.itau_EmitterRecombination(), '--', label='Emitter')

    def PlotComponants(self, plt=False):
        if plt == False:
            plot(self.Deltan, self.tau_SRH(), '--', label='SRH')
            plot(self.Deltan, self.tau_aug(), '--', label='Auger')
            plot(self.Deltan, self.tau_rad(), '--', label='Radiative')
            plot(self.Deltan, self.tau_EmitterRecombination(),
                 '--', label='Emitter')
        else:
            plt.plot(self.Deltan, self.tau_SRH(), '--', label='SRH')
            plt.plot(self.Deltan, self.tau_aug(), '--', label='Auger')
            plt.plot(self.Deltan, self.tau_rad(), '--', label='Radiative')
            plt.plot(
                self.Deltan, self.tau_EmitterRecombination(), '--', label='Emitter')

    def iTotalLifetime(self):

        return (self.itau_aug() + self.itau_EmitterRecombination() + self.itau_rad() + self.itau_SRH())

    def TotalLifetime(self):

        return 1. / self.iTotalLifetime()

    def Plot(self, plt=False):
        if plt == False:
            plt = figure('Tau')
        ax = plt.add_subplot(111)

        ax.plot(self.Deltan, self.TotalLifetime(), 'k-', label='Lifetime')

        legend(loc=0)
        loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('Lifetime (s)')

    def iPlot(self, plt=False):
        if plt == False:
            plt = figure('i-Tau')
        ax = plt.add_subplot(111)

        ax.plot(self.Deltan, self.iTotalLifetime(), '-', label='Lifetime')

        legend(loc=0)
        loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('Inverse Lifetime (s$^{-1}$)')


class fit_tau(SRH, Auger, Radiative):

    def __init__(self, tau_m, Deltan):
# Constants.__init__(self)       #not needed?, as inherrants properties
# from SHR,Ayger and Radiative it gets the Constants Class
        SRH.__init__(self)
        Auger.__init__(self)
        Radiative.__init__(self)
        self.tau_m = tau_m
        self.Deltan = Deltan  # overwrites default value in class "Constants"

    def Worst_SRV_DoubleSide(self):
        # This assumes the bulk lifetime is Auger limited

        self.Deltan = self.Deltan[~isnan(self.tau_m)]
        self.tau_m = self.tau_m[~isnan(self.tau_m)]

        S = self.Width * (1 / self.tau_m - self.itau_aug()) / 2

        return S

    def fitting(self):
        plt = figure()
        ax1 = plt.add_subplot(121)
        ax2 = plt.add_subplot(122)

        # cleaning up input signal
        self.Deltan = self.Deltan[~isnan(self.tau_m)]
        self.tau_m = self.tau_m[~isnan(self.tau_m)]

# itau = 1/self.tau - 1/self.tau_aug(Deltan)     #xxx fix me
# tau_SRH = 1/itau_SRH    #tau_SRH, ambiguous name, irection

        # limits fit to exclude Et, why is this done? why not just reference
        # tau_SRH?
        f1 = lambda Deltan, tau_n, tau_p: self.tau_SRH(
            Deltan, abs(tau_n), abs(tau_p))
        f2 = lambda Deltan, tau_n, tau_p: self.itau_SRH(
            Deltan, abs(tau_n), abs(tau_p))

        tau_SHR = (1 / self.tau_m - self.itau_aug() - self.itau_rad())**(-1)
        itau_SHR = 1 / self.tau_m - self.itau_aug() - self.itau_rad()

        popt, pcov = scipy.optimize.curve_fit(
            f1, self.Deltan, tau_SHR, p0=[1e-3, 1e-3])
        ipopt, ipcov = scipy.optimize.curve_fit(
            f2, self.Deltan, itau_SHR, p0=[1e-3, 1e-3])

        # dummy Et vals
        popt = numpy.append(popt, 0)
        ipopt = numpy.append(ipopt, 0)
#        print 'popt', popt
#        print 'ipopt', ipopt

        ax1.plot(self.Deltan, self.tau_m, label='Original')
#        ax1.plot(self.Deltan,self.tau_SRH(Deltan),'-',label='tau srh')
        ax1.plot(self.Deltan, 1 / (self.itau_SRH(self.Deltan, *popt) +
                                   self.itau_aug() + self.itau_rad()), '--', label='Fitted')

        ax2.plot(self.Deltan, 1 / self.tau_m / 1e3, label='Original')
        ax2.plot(self.Deltan, (self.itau_SRH(self.Deltan, *ipopt) +
                               self.itau_aug() + self.itau_rad()) / 1e3, '--', label='Fitted')
#        ax2.plot(self.Deltan,itau_SRH,label='itau_SRH')

        print '\t\t\t  tau_n \t\t tau_p \t\t E_t'
        print 'Inverse tau fitting \t {0:.2e} \t {1:.2e}\t {2:.2e}'.format(ipopt[0], ipopt[1], ipopt[2])
        print 'tau fitting \t\t {0:.2e} \t {1:.2e}\t {2:.2e}'.format(popt[0], popt[1], popt[2])

        ax1.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax2.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax1.set_ylabel('Lifetime (s)')
        ax2.set_ylabel('Inverse Lifetime (ms$^{-1}$)')

        ax1.loglog()
        ax2.loglog
        legend(loc=0)
        # show()


def CheckFittingJ0():

    # Constants.Deltan=logspace(12,20,1000)
    # Constants.Deltan=logspace(12,20,1000)

    #Deltan = Constants.Deltan
    #Deltan = Constants.Deltan
    Constants.Deltan = logspace(14, 17, 1000)
    plt2 = figure()
    plt = figure()

    # Auger().PlotAll(plt)
    # Radiative().PlotAll(plt)
    # SRH().PlotAll(plt)

    TotalLifetime().J0e = 3e-13

    TotalLifetime().iPlot(plt2)

    # Emitter_Recombination.J0e = 1e-13
    # Emitter_Recombination().PlotAll(plt)
    # TotalLifetime().Plot(plt)
    text(1e12, 1e-1, 'Doping: ' + str(Constants.Doping) +
         ' cm^-3', ha='center', va='center')
    text(1e12, 0.7e-1, 'Width: ' + str(Constants.Width) +
         ' cm', ha='center', va='center')
    text(1e12, 0.4e-1, 'J0e: ' + str(Emitter_Recombination().J0e) +
         ' Acm^-2', ha='center', va='center')

    yscale('linear')

    TotalLifetime().iPlot(plt2)
    yscale('linear')
    xscale('linear')
    show()


#A = Auger()
#R = Radiative()
#A.Doping= 8e15
#R.Doping = 8e15
# A.T = 300
# R.T = 300

# R.tau_rad()
# print R.Doping


# plot(A.Deltan,1./(A.itau_aug()+R.itau_rad()))
# savetxt('Intrinsic',vstack((A.Deltan,1/(A.itau_aug()+R.itau_rad()))).T,delimiter='\t')
# plot(R.Deltan,R.tau_rad(),'--')
# plot(R.Deltan,A.tau_aug(),'-.')

# print A.tau_aug()[0]*1e6,1./(A.itau_aug()+R.itau_rad())[0]
# print R.tau_rad()[0]
# loglog()


# semilogx()
# figure()
# plot(A.Deltan,A.itau_aug()+R.itau_rad())
# plot(R.Deltan,R.itau_rad(),'--')
# plot(R.Deltan,A.itau_aug(),'-.')
# loglog()
# show()
''''
data=genfromtxt('.data_for_RecombinationFile',skip_header=2)
index = where(data[:,2]==max(data[:,2]))[0]
data = data[index:]
index = where( data[:,0] > 0)
data = data[index]

Constants.Doping=1e16

#a =fit_tau()

wafer1 = fit_tau(data[:,0],data[:,2])
wafer1.Width= 0.035
wafer1.fitting()

figure()

for i in linspace(1e16,5e16,5):
    This = Constants()
    Constants.Doping = i
    print Constants.Doping
    print Auger().Doping
    itau = Auger().itau_aug()+SRH().itau_SRH(Constants.Deltan,1e5,0,0)+Radiative().itau_rad()


    plot(Constants.Deltan,1/itau,label = str(i))
    plot(data[:,2],data[:,0],'.')


legend(loc=0)
loglog()

Constants.Doping=1e16

figure('SVR')
plot(wafer1.Deltan,wafer1.Worst_SRV_DoubleSide())
loglog()

show()
'''
'''
def GetDeltanGivenTau(W):

    return  a.Deltan[abs(a.Deltan/a.tau_aug()-G/W*ones(a.Deltan.shape[0])).argmin(axis=0)]

a = Auger()
G = 2.5e17
Width = 0.001
wide =array([])
carrierdensity =array([])
#print store
for Width in logspace(-4,-1,100):
    wide = append(wide, array(Width))
    carrierdensity = append(carrierdensity,array(GetDeltanGivenTau(Width)))
    #a.Width = W
    #a.Deltan = 1e16
    #for i in range(50):

        #tau = a.tau_aug()
        #a.Deltan = G/a.Width*tau
        #print tau,a.Width
    #Voltage  = 0.02585 *2*log(G/w *
#print a.Deltan[abs(a.Deltan/a.tau_aug()-G/Width*ones(a.Deltan.shape[0])).argmin(axis=0)]
#plot(a.Deltan,abs(a.Deltan/a.tau_aug()-G/Width*ones(a.Deltan.shape[0])))
#plot(a.Deltan,a.Deltan/a.tau_aug())
#plot(a.Deltan,G/0.001*ones(a.Deltan.shape[0]))
print wide.shape
plot(wide,2*0.02585*log(carrierdensity/a.ni),'.')
#plot(wide,2*0.02585*log(carrierdensity/a.ni/a.ni))
semilogx()
#semilogy()
'''
# show()

"""
FileName = 'F-FZL.txt'
Dn = genfromtxt(FileName,usecols=(8,),delimiter='\t',skip_header=18,skip_footer=29)
tau_m = genfromtxt(FileName,usecols=(11,),delimiter='\t',skip_header=18,skip_footer=29)

wafer1 = fit_tau(tau_m,Dn)
wafer1.fitting()
"""
