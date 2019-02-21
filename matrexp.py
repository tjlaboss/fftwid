import scipy as sp
import scipy.linalg as la
import depletr
from depletr.nuclides.thermal import u235, u238, pu239

test_matrix = sp.array([
# _from_   U235             U238            Pu239        |_to_
	[ u235.sigma_a,               0,               0],  # U235
	[              0,  u238.sigma_a,               0],  # U238
	[              0, -u238.sigma_y,   pu239.sigma_a]   # Pu239
], dtype=float)


NSTEPS = 10
num = sp.zeros((3, NSTEPS))
num[:, 0] = (0.0495, 0.95, 0.0005)
flux = 1
dt = 0.1
for k in range(NSTEPS):
	t = dt*k
	dn = la.expm(-test_matrix*flux*dt)
	if k > 0:
		num[:, k] = num[:, k-1].dot(dn)
	print('time t:', t, '\n', num[:, k])
	print('dn\n', dn, '\n')
	
	
PLOT = True
if PLOT:
	from pylab import *
	pf = semilogy
	#pf = plot
	pf(num[0, :], label=u235.name)
	pf(num[1, :], label=u238.name)
	pf(num[2, :], label=pu239.name)
	legend()
	show()


