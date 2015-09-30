

import recombination.SRH as test
import matplotlib.pylab as plt
import numpy as np 


fig, ax = plt.subplots(1)
a = test.SRH()
defect = test.SHR_defects(a)


a.Deltan = np.logspace(12,17)
print a.tau_n, '\n\n\n'

defect.Fe(1e10)
a.PlotAll(fig)

defect.FeB(1e10)
a.PlotAll(fig)

plt.show()

print 'hello'
