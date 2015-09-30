import numpy
from pylab import *
import BlackBody
from glob import glob
import sys
import os
import scipy.constants as Const


class IntrinsicCarrierDensity():
    # This is taken from http://www.azom.com/article.aspx?ArticleID=8489#4
    ni = 18.9e-18  # cm-3