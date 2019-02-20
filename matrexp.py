import scipy as sp
import scipy.linalg as la


from depletr.data import u235_t, u238_t, pu239_t

test_matrix = sp.array([
# _from_   U235             U238            Pu239         | _to_
	[ u235_t.sigma_a,               0,               0],  # U235
	[              0,  u238_t.sigma_a,               0],  # U238
	[              0, -u238_t.sigma_y, pu239_t.sigma_a]   # Pu239
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
	pf(num[0, :], label=u235_t.name)
	pf(num[1, :], label=u238_t.name)
	pf(num[2, :], label=pu239_t.name)
	legend()
	show()


