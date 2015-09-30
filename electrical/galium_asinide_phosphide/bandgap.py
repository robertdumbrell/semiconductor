import numpy
from pylab import *
import BlackBody
from glob import glob
import sys
import os
import scipy.constants as Const


class Bandgap():
    # This is taken from http://www.azom.com/article.aspx?ArticleID=8489#4
    E = 1.98 * Const.e  # eV