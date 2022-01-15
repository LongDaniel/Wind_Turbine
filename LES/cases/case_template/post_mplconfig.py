import numpy as np
import matplotlib as mpl

# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

fig_width_pt = 550.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
  'axes.labelsize': 12,
  'text.fontsize': 10,
  'legend.fontsize': 10,
  'xtick.labelsize': 10,
  'ytick.labelsize': 10,
  'text.usetex': True,
  'figure.figsize': fig_size}
mpl.rcParams.update(params)
mpl.rcParams['text.latex.preamble'].append(r'\usepackage{amsmath}')

