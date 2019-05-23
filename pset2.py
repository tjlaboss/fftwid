# 22.215, PSet02
#   by Travis Labossiere-Hickman (email: tjlaboss@mit.edu)
#
# Questions 1-4

import depletr
import scipy as sp


SPECTRUM = "fast"
AVERAGE_KINF = {"thermal": 1.3169320996513936,
				"fast"   : 0.9892984136969714}
ELEMENTS = ("Pu", "Am", "Cm")
# found with quicksolve.py
MOX_AT = sp.array([0, 0.03182388222963699, 0.03319650338671934, 0.03356631994343097])
DESCRIPTOR = ("Fresh      ", "Once-Burnt ", "Twice-Burnt", "Thrice-Burnt")
NSTEPS = 100

study = depletr.Depleter(power=35, enrichment=5, max_burnup=50, spectrum=SPECTRUM)
all_nuclides = study.get_all_nuclides()
num_nuclides = len(all_nuclides)
cycle_wts = sp.zeros((num_nuclides + 2, len(MOX_AT)))

cout = None
for i, mox_frac in enumerate(MOX_AT):
	print(DESCRIPTOR[i], "Fuel")
	print("Mox fraction: {:6.3%}".format(mox_frac))
	if i == 0:
		# Fresh fuel
		cout = study.deplete_fresh(NSTEPS, plots=0)
	else:
		cin = study.reprocess(cout, ELEMENTS, mox_frac)
		print(100*cin/cin.sum())
		cout = study.reload(cin, NSTEPS, verbose=True)
	cycle_wts[:, i] = depletr.fuel.at_to_wt_arbitrary(all_nuclides, cout)
	print()


fmt_dict = {'float_kind':"{:11.5e}".format}
print('\nWeight Fraction at EOC\n' + " "*32, end=" ")
for desc in DESCRIPTOR:
	print(desc, end=" "*2)
print()
for j, nuc in enumerate(all_nuclides):
	cycstr = sp.array2string(cycle_wts[j], formatter=fmt_dict, separator="  ")
	nucstr = str(nuc)
	print(nucstr + " "*(32-len(nucstr)) + cycstr)
