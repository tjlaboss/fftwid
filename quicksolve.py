import depletr
import scipy as sp
from scipy.optimize import fsolve

SPECTRUM = "fast"
AVERAGE_KINF = {"thermal": 1.3169320996513936,
				"fast"   : 0.9892984136969714}
ELEMENTS = ("Pu", "Am", "Cm")
MOX_AT = {"fast"   : sp.array([0, 0.03182388222963699, 0.03319650338671934, 0.03356631994343097]),
          "thermal": sp.array([0, 0.2987412473925805])}
DESCRIPTOR = ("Fresh", "Once-Burnt", "Twice-Burnt", "Thrice-Burnt")
NSTEPS = 100

study = depletr.Depleter(power=35, enrichment=5, max_burnup=50, spectrum=SPECTRUM)
cout = study.deplete_fresh(NSTEPS, plots=0, verbose=False)


# Once-burnt fuel
cin = study.reprocess(cout, ELEMENTS, mox_frac=MOX_AT[SPECTRUM][1])
cout = study.reload(cin, NSTEPS)

# Twice-burnt fuel
cin = study.reprocess(cout, ELEMENTS, mox_frac=MOX_AT[SPECTRUM][2])
cout = study.reload(cin, NSTEPS)
# Thrice-burnt fuel
#cin = study.reprocess(cout, which_elements=ELEMENTS, mox_frac=MOX_AT[SPECTRUM][3])

def optimize_mox_mix(n):
	global cout
	cin = study.reprocess(cout, which_elements=ELEMENTS, mox_frac=n)
	kinf3 = study.reload(cin, NSTEPS, return_kinf=True)
	print(n)
	return kinf3.mean()


n3 = fsolve(lambda n: optimize_mox_mix(n) - AVERAGE_KINF[SPECTRUM], 0.01)[0]
print("n3:", n3)
print("kinf:", optimize_mox_mix(n3))

