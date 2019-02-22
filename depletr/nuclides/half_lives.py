# Half-Lives
#
# Source: NNDC.BNL.gov/nudat2

from scipy import log

LN2 = log(2)
MINUTE = 60
HOUR = 3600
DAY = 24*HOUR
YEAR = 365.25*DAY

# Alpha decay
ALPHA = {
	# Uranium
	"U235"  : 703.8E6*YEAR,
	"U236"  : 2.342E7*YEAR,
	# Neptunium
	"Np237" : 2.144E6*YEAR,
	# Plutonium
	"Pu238" : 87.7*YEAR,
	"Pu239" : 24110*YEAR,
	"Pu240" : 6561*YEAR,
	"Pu241" : 14.29*YEAR,
	"Pu242" : 3.73E5*YEAR,
	# Americium
	"Am241" : 432.6*YEAR,
	"Am242m": 141*YEAR/0.0046,
	"Am243" : 7364*YEAR,
	# Curium
	"Cm242" : 162.86*DAY,
	"Cm243" : 28.9*YEAR,  # but 99.71%
	"Cm244" : 18.11*YEAR,
	"Cm245" : 8453*YEAR
}

# Beta-
BETAM = {
	# Thorium
	"Th233" : 22.3*MINUTE,
	"Th244" : 24.1*DAY,
	# Palladium
	"Pa232" : 1.31*DAY,
	"Pa233" : 26.975*DAY,
	"Pa234" : 1.159*MINUTE/.9984,
	# Uranium
	"U237"  : 6.75*DAY,
	"U239"  : 23.45*MINUTE,
	# Neptunium
	"Np236" : 22.5*HOUR/.5,
	"Np238" : 2.117*DAY,
	"Np239" : 2.365*DAY,
	"Np240" : 7.22*MINUTE,
	# Americium
	"Am242" : 16.02*HOUR/.872,
}

# Beta+ (and electron capture)
BETAP = {
	# Neptunium
	"Np236" : 22.5*HOUR/.5,
	# Americium
	"Am242" : 16.02*HOUR/.173,
}

# Gamma (only for Americium metastable)
GAMMA = {"Am242m": 141*YEAR/0.9954}
