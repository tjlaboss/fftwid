# 22.215, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Question 5

import depletr
from scipy.constants import N_A

AVG_MASS = 239
MEV_PER_FISSION = 200
MEV_PER_MWD = 5.39266E23
NHM = 1000*N_A/AVG_MASS
X = MEV_PER_FISSION/MEV_PER_MWD
BU = 0.15*NHM*X

ELEMENTS = ("Pu", "Am", "Cm")
NSTEPS = 500

# PWR with fresh fuel
thermal_deplet = depletr.Depleter(power=35, enrichment=5, max_burnup=50, spectrum="thermal")
cout = thermal_deplet.deplete_fresh(NSTEPS, plots=0, verbose=False)
# Fast reactor with once-burnt fuel
cin = thermal_deplet.reprocess(cout, ELEMENTS, mox_frac=0.3)
fast_deplet = depletr.Depleter(power=35, enrichment=-0.71, max_burnup=BU, spectrum="fast")
cout = fast_deplet.reload(cin, NSTEPS)

print(cout)
