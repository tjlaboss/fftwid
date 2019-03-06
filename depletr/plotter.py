# Plotter

import numpy as np
import matplotlib.pyplot as plt
from depletr.fuel import wt_to_at_uranium
from nuclides.half_lives import MINUTE, HOUR, DAY, MONTH, YEAR

UCOLOR = "darkgoldenrod"  # because of the name
FLUXCOLOR = (0.35, 0.35, 0.35)  # Dark gray
UNITS = {"second": 1,
         "minute": MINUTE,
         "hour"  : HOUR,
         "day"   : DAY,
         "month" : MONTH,
         "year"  : YEAR}


def _get_axis(ax):
	if ax is None:
		ax = plt.figure().add_subplot()
	return ax


def _get_t_in_unit(tvals, unit):
	basestr = "$t$ ({}s)"
	if unit not in UNITS:
		print("Invalid unit", unit)
		return tvals, basestr.format("second")
	return tvals/UNITS[unit], basestr.format(unit)


def make_actinides_plot(tvals, num, all_nuclides, ax=None, unit="day",
                        plot_type="semilogy", deadend_actinides=False, fission_products=False):
	ax = _get_axis(ax)
	tvals, tstr = _get_t_in_unit(tvals, unit)
	plot_f = {"semilogy": ax.semilogy,
	          "semilogx": ax.semilogy,
	          "loglog"  : ax.loglog,
	          "plot"    : ax.plot}[plot_type]
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
	ax.set_xlabel(tstr)
	ax.grid(True, which="both", ls="-")
	ax.set_xlim(tvals[0], tvals[-1])
	ax.legend(fancybox=True, shadow=True, bbox_to_anchor=(-0.13, 1.0))
	ax.set_title("Actinide Number Density", fontweight="bold")
	return ax


def make_enrichment_flux_plot(tvals, enrichvals, fluxvals, ax=None, unit="day"):
	ax = _get_axis(ax)
	tvals, tstr = _get_t_in_unit(tvals, unit)
	# Enrichment
	handles = []
	labels = []
	wtfrac = np.zeros(enrichvals.shape)
	for i, e in enumerate(enrichvals):
		wtfrac[i] = wt_to_at_uranium(e)
	ax.plot(tvals, wtfrac*100, color=UCOLOR, label="$^{235}$U")
	hu, lu = ax.get_legend_handles_labels()
	handles += hu
	labels += lu
	ax.yaxis.label.set_color(UCOLOR)
	ax.spines['right'].set_color(UCOLOR)
	ax.set_ylabel("Uranium 235 Enrichment (wt%)")
	# Flux
	axf = ax.twinx()
	axf.plot(tvals, fluxvals*1E24, color=FLUXCOLOR, linewidth=2, label="$\phi$")
	hf, lf = axf.get_legend_handles_labels()
	handles += hf
	labels += lf
	axf.set_xlabel(tstr)
	axf.yaxis.label.set_color(FLUXCOLOR)
	[t.set_color(FLUXCOLOR) for t in axf.yaxis.get_ticklines()]
	[t.set_color(FLUXCOLOR) for t in axf.yaxis.get_ticklabels()]
	axf.spines['right'].set_color(FLUXCOLOR)
	axf.set_ylabel("$\phi(t)$ (neutrons/cm${}^2$/s)")
	ax.grid(True, which="both", ls="-")
	ax.legend(handles, labels, loc="center left")
	ax.set_title("Enrichment and Flux", fontweight="bold")
	return ax


def make_element_depletion_plot(tvals, list_of_arrays, list_of_nuclides, ax=None, unit="day",
                                relative=True, element=None):
	ax = _get_axis(ax)
	tvals, tstr = _get_t_in_unit(tvals, unit)
	if len(list_of_arrays) != len(list_of_nuclides):
		print("Mismatch in length of `list_of_arrays` and `list_of_nuclides`.")
		return ax

	for i, nuclide in enumerate(list_of_nuclides):
		nucvals = list_of_arrays[i]
		if relative:
			nucvals /= nucvals[0]
		color = UCOLOR if i == 0 else None
		ax.plot(tvals, nucvals, color=color, label=nuclide.latex)
	ax.grid()
	ax.set_xlabel(tstr)
	titstr = ""
	if relative:
		titstr += "Relative "
	if element is None:
		element = list_of_nuclides[0].element
	titstr += element
	titstr += " Concentration"
	ax.set_title(titstr, fontweight="bold")
	ax.legend()
	return ax


def make_kinf_plot(tvals, kvals, ax=None, unit="day"):
	ax = _get_axis(ax)
	tvals, tstr = _get_t_in_unit(tvals, unit)
	ax.plot(tvals, kvals, 'ro', label="$k_\infty$")
	ax.set_xlabel(tstr)
	#ax.set_ylabel("$k_\infty$", fontsize=13, color='r')
	if kvals.min() > 1:
		ax.set_ylim(1, kvals.max() + 0.1)
	ax.grid()
	ytickvals = np.arange(np.floor(10*kvals.min()), np.ceil(10*kvals.max()+1), 1)/10
	ax.set_yticks(ytickvals)
	ax.set_title("k-infinity", fontweight="bold")
	ax.legend()
	return ax
