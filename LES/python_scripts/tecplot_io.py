import numpy as np

def tecplot_reader(file, ndim, nheaderline):
  """Tecplot reader."""
  arrays = []
  with open(file, 'r') as a:
    for idx, line in enumerate(a.readlines()):
      if idx < nheaderline:
        continue
      else:
        arrays.append([float(s) for s in line.split()])

  arrays = np.concatenate(arrays)
  output = arrays.reshape(ndim)
  #output = arrays.copy()
  #output = np.split(arrays, nb_var)

  #for var in output:
  #  print(var.shape)
  #  var.reshape(ndim)

  return output
