# Printer

import numpy as np


def print_decay_results(names, times, concentrations, other_actinides="other"):
	assert concentrations.shape == (len(names)+1, len(times)+1), \
		"Misshapen data. Format is: concentrations[nuclides+1, times+1]"
	names += [other_actinides]
	ndec = 16
	nchar = len(max(names, key=len))
	print("\n\nDecay after __ years:\n")
	
	fmt_header = '{:^' + str(ndec) + '}|'
	headstr = " "*nchar + " | "
	for t in [0] + list(times):
		headstr += fmt_header.format(t)
	print(headstr)
	print("-"*len(headstr))
	datastr = np.array2string(concentrations, max_line_width=np.inf, separator=' | ')
	datastr = datastr.replace('[', ' ').replace(']', ' ').split('\n')
	fmt_nuclide = '{:^' + str(nchar) + '} |'
	for i, name in enumerate(names):
		print(fmt_nuclide.format(name) + datastr[i])
	print()
