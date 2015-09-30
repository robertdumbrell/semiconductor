import numpy
import scipy.constants as Const
import os



class OpticalProperties_GaAsP():

    alpha_version = 'Brianna2014'
    n_version = ''

    def __init__(self):
        self.initalise_optical_constants()

    def initalise_optical_constants(self):

        if os.path.isdir("C:\Users\mattias"):
            self.Folder = r'C:\Users\mattias'
        elif os.path.isdir("C:\Users\z3186867"):
            self.Folder = r'C:\Users\z3186867'
        else:
            print 'You need to find a new location for the use directoy'

        self.Alpha_Folder = r'\Dropbox\CommonCode\SpectralPL\Matterials Absorption Data_PVlighthouse'
        self.load_optical_properties()
