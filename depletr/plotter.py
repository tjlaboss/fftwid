# Plotter

import numpy as np
import matplotlib.pyplot as plt
from depletr.fuel import wt_to_at_uranium

UCOLOR = "darkgoldenrod"  # because of the name
FLUXCOLOR = (0.35, 0.35, 0.35)  # Dark gray


def _get_axis(ax):
	if ax is None:
		ax = plt.figure().add_subplot()
	return ax


def make_actinides_plot(tvals, num, all_nuclides, ax=None, plot_f=plt.semilogy,
						deadend_actinides=False, fission_products=False):
	ax = _get_axis(ax)
	# Actinides
	for i, nuclide in enumerate(all_nuclides):
		if i == 0:
			plot_f(tvals, num[i], color=UCOLOR, linewidth=2, label=nuclide.latex)
		else:
			plot_f(tvals, num[i], label=nuclide.latex)
	if deadend_actinides:
		plot_f(tvals, num[-2], ':', label="other")
	if fission_products:
		plot_f(tvals, num[-1], '-', label="PRODUCT")
	ax.set_xlabel("$t$ (days)")
	ax.grid(True, which="both", ls="-")
	ax.set_xlim(tvals[0], tvals[-1])
	ax.legend(fancybox=True, shadow=True, bbox_to_anchor=(-0.13, 1.0))
	ax.set_title("Actinide Number Density", fontweight="bold")
	return ax


def make_enrichment_flux_plot(tvals, enrichvals, fluxvals, ax=None):
	ax = _get_axis(ax)
	# Enrichment
	handles = []
	labels = []
	wtfrac = np.zeros(enrichvals.shape)
	for i, e in enumerate(enrichvals):
		wtfrac[i] = wt_to_at_uranium(e)
	plt.plot(tvals, wtfrac*100, color=UCOLOR, label="$^{235}$U")
	hu, lu = ax.get_legend_handles_labels()
	handles += hu
	labels += lu
	ax.yaxis.label.set_color(UCOLOR)
	ax.spines['right'].set_color(UCOLOR)
	ax.set_ylabel("Uranium 235 Enrichment (wt%)")
	# Flux
	axf = plt.twinx()
	axf.plot(tvals, fluxvals*1E24, color=FLUXCOLOR, linewidth=2, label="$\phi$")
	hf, lf = axf.get_legend_handles_labels()
	handles += hf
	labels += lf
	axf.set_xlabel("$t$ (days)")
	axf.yaxis.label.set_color(FLUXCOLOR)
	[t.set_color(FLUXCOLOR) for t in axf.yaxis.get_ticklines()]
	[t.set_color(FLUXCOLOR) for t in axf.yaxis.get_ticklabels()]
	axf.spines['right'].set_color(FLUXCOLOR)
	axf.set_ylabel("$\phi(t)$ (neutrons/cm${}^2$/s)")
	ax.grid(True, which="both", ls="-")
	plt.legend(handles, labels, loc="center left")
	plt.title("Enrichment and Flux", fontweight="bold")


def make_element_depletion_plot(tvals, list_of_arrays, list_of_nuclides, ax=None, relative=True, element=None):
	ax = _get_axis(ax)
	if len(list_of_arrays) != len(list_of_nuclides):
		print("Mismatch in length of `list_of_arrays` and `list_of_nuclides`.")
		return ax

	for i, nuclide in enumerate(list_of_nuclides):
		nucvals = list_of_arrays[i]
		if relative:
			nucvals /= nucvals[0]
		color = UCOLOR if i == 0 else None
		plt.plot(tvals, nucvals, color=color, label=nuclide.latex)
	plt.grid()
	titstr = ""
	if relative:
		titstr += "Relative "
	if element is None:
		element = list_of_nuclides[0].element
	titstr += element
	titstr += " Concentration"
	plt.title(titstr, fontweight="bold")
	plt.legend()

