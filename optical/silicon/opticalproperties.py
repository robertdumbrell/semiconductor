
import numpy as np
import sys
import os
import scipy.constants as Const


class OpticalProperties():

    """
    A class containg the optical constants of silicon
    These are temperature dependence.
    """

    def __init__(self,  matterial='Si', alpha_version='Schinke2014', n_version='Green2008'):
        print 'Does not currently take matterial inputs'
        self.alpha_version = alpha_version
        self.n_version = n_version
        self.initalise_optical_constants()

    # something else
    def initalise_optical_constants(self):
        """
        Sets the folder to find tabularted data and loads it
        """

        self.Folder = os.path.dirname(os.path.realpath(__file__))
        self.load_optical_properties()

    def wavelength_to_alpha(self, wavelength):
        return np.interp(wavelength,
                         self.wavelength_emission,
                         self.alpha_BB)

    def wavelength_to_refactiveindex(self, wavelength):
        return np.interp(wavelength,
                         self.wavelength_emission,
                         self.n)

    def load_optical_properties(self, alpha_version=False, n_version=False):
        """
        Loads alpha and n
        from the provided or from self. the name
        """

        if not alpha_version:
            alpha_version = self.alpha_version
        else:
            self.alpha_version = alpha_version
        # print 'here', self.alpha_version
        if not n_version:
            n_version = self.n_version
        else:
            self.n_version = n_version

        # Getting Alpha
        data = getattr(self, 'optical_' + alpha_version)()
        self.wavelength_emission, self.alpha_BB = data[
            'Wavelength'], data['alpha']

        # GEtting n
        data = getattr(self, 'optical_' + n_version)()
        # adjust n to the same wavelength as, what ever alpha is at, and then
        # save it
        self.n = np.interp(
            self.wavelength_emission, data['Wavelength'], data['n'])
        self.k = np.interp(
            self.wavelength_emission, data['Wavelength'], data['k'])

    def change_temp_Green2008(self, Temp=False):
        """
        Uses green2008 data to account for temp changes in optical constants
        Uses a model he proped   alpha propt Temp^coef
        """

        if not Temp:
            Temp = self.Temp

        self.load_optical_properties()

        Folder = os.path.dirname(os.path.realpath(__file__))
        File = r'\Silicon_Green08'
        data = genfromtxt(
            Folder + File, names=True, delimiter='\t', filling_values=0)

        wavelength = data['Wavelength'] * 1000
        tempcof = data['C_ka']
        # Adjusting to the emission values
        tempcof = np.interp(self.wavelength_emission, wavelength, tempcof)
        # Then we correct the alpha values for temp
        b = tempcof * 300 * 1e-4  # *self.alpha_BB
        self.alpha_BB *= (Temp / 300.)**b
        self.k *= (Temp / 300.)**b

        # Now we correct n, C_n
        tempcof = data['C_n']
        # Adjusting to the wavelegnth range used
        tempcof = np.interp(self.wavelength_emission, wavelength, tempcof)
        # Finding the temp cof
        b = tempcof * 300 * 1e-4  # *self.alpha_BB
        self.n *= (Temp / 300.)**b

    def optical_Daub1995(self):
        """
        wavelegnth and alpha from photoluminesence
        """
        # Only has alpha
        self.Alpha_FileName = '\Silicon_Daub'
        data = np.genfromtxt(self.Folder
                             + self.Alpha_FileName, names=True)
        return data

    def optical_Green2008(self):
        """
        combined sets of alpha, n and temperature coeffiecnts. Uses Konges krange to correct k
        """
        # HAs alpha, n, k
        self.Alpha_FileName = '\Silicon_Green08'
        data = np.genfromtxt(self.Folder
                             + self.Alpha_FileName, delimiter="\t", names=True)
        data['Wavelength'] *= 1000
        return data

    def optical_Green1995(self):
        """
        combined sets of alpha, n and temperature coeffiecnts
        """
        # HAs alpha, n, k
        self.Alpha_FileName = '\Silicon_Green95'
        data = np.genfromtxt(self.Folder
                             + self.Alpha_FileName, delimiter="\t", names=True)
        return data

    def optical_Schinke2014(self):
        """
        wavelength and alpha
        Combination of well chacterised reflection and transmission and photoluminescence measurements
        Comes with uncertainty values on each reported value
        Planning follow up paper with comparson to EQE measurements.
        """
        # Only has alpha
        self.Alpha_FileName = '\Silicon_Schinke2014'
        data = np.genfromtxt(self.Folder
                             + self.Alpha_FileName, delimiter="\t", names=True)
        return data

    def optical_Schinke2015(self):
        """
        wavelength and alpha
        Combination of well chacterised reflection and transmission and photoluminescence measurements
        Combines  R&T, ellipometry, EQE and PL measurements at 295K
        """
        # Only has alpha
        self.Alpha_FileName = '\Silicon_Schinke2015'
        data = np.genfromtxt(self.Folder
                             + self.Alpha_FileName, delimiter="\t", names=True)
        return data
