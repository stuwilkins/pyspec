from pylab import *
from pyspec import spec, fit, fitfuncs

def fittest():
    sf = spec.SpecDataFile('/Users/swilkins/Dropbox/Data/X1A2/oct09_2/Mn214_Surface/Mn214_Surface.01')
    scan = sf[229]
    f = fit.fit(x = scan.H, y = scan.Detector, funcs = [fitfuncs.linear, fitfuncs.gauss],
                optimizer = 'levmar')
    f.go()

if __name__ == '__main__':
    fittest()
