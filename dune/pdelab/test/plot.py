from pandas import *
import matplotlib.pyplot as plt
import itertools
import os
import subprocess
import re
from collections import OrderedDict
#from slurmpy import Slurm

def tuple_to_string(tpl):
  return "_".join([str(x) for x in tpl])

def run (bincall, plotprefix, configuration, permutations, extractions, xlabel = "", ylogscale = False):

  results = {}
  for tpl in list(itertools.product(*permutations.values())):
    results[tpl[1:]] = {}
    for extraction in extractions.keys():
      results[tpl[1:]][extraction] = {}


  for tpl in list(itertools.product(*permutations.values())):

    runconfig = configuration.copy()
    runconfig.update(zip(permutations.keys(), tpl))
    params = ' '.join(['-%s %s' % (key, value) for (key, value) in runconfig.items()])

    print ("Running " + str(tpl) + "\n\t" + params)

    #sbatch --out "test.out" batch.sbatch "params..."
    use_slurm = False
    if use_slurm:
        os.system("sbatch --job-name=prmstudy --out \"out/" + tuple_to_string(tpl) + ".out\" batch.sbatch \"" + params + "\"")
    else:
        os.system("mpirun -np 20 ./testgeneo " + params + " > 'out/" + tuple_to_string(tpl)+".out'")


    #os.system (bincall + " " + params + " > out")

  if use_slurm:
      print("waiting for runs...")
      os.system("while squeue | grep --quiet prmstudy; do sleep 3; done")

  for tpl in list(itertools.product(*permutations.values())):
    print ("Reading out/" + tuple_to_string(tpl) + ".out")

    with open("out/" + tuple_to_string(tpl) + ".out", encoding="utf-8") as f:
      for line in f:
        #print(line)
        for extraction in extractions.keys():
          match = re.search('(?<=' + extractions[extraction] + ').*', line)
          if (match):
            results[tpl[1:]][extraction][tpl[0]] = match.group(0)



  primaries = list(permutations.values())[0]

  for extraction in extractions.keys():

    #fig, ax = plt.plot()

    fig = plt.figure()

    print (permutations.values())
    secondaries = list(permutations.values())
    secondaries.pop(0)
    for subtpl in list(itertools.product(*secondaries)):
      line, = plt.plot(np.array(primaries), [float(results[subtpl][extraction].get(x, 0.0)) for x in primaries])
      line.set_label(" ".join(subtpl))

    #ax.legend()
    plt.ylabel(extraction)
    if (xlabel == ""):
      plt.xlabel(list(permutations.keys())[0])
    else:
      plt.xlabel(xlabel)
    #plt.title(extraction)


    if ylogscale:
      plt.yscale('symlog')
    #plt.xscale('symlog')


    plt.legend()
    fig.savefig(plotprefix + extraction + ".svg")
    fig.clf()

    plt.close()
  #plt.show()











config = {
  "cells": "250",
  "subdomainsx": "10",
  "subdomainsy": "1",
  "nev": "20",
  "nev_arpack": "20",
  "overlap": "1",
  "layers": "25",
  "contrast": "1e5",
  "method": "geneo",
  "part_unity": "standard",
  "hybrid": "false",
  "coarse_only": "false",
  "arpack_factor": "1.0",
  "dirichlet_interpol_in_coarse_only": "false",
}

config_approx = config.copy()
config_approx["coarse_only"] = "true"
config_approx["subdomainsx"] = "4"
config_approx["subdomainsy"] = "4"

config["layers"] = "40" # FIXME
config["cells"] = "400" # FIXME
config["subdomainsx"] = "10" # FIXME
config["subdomainsy"] = "2"
#config_approx["layers"] = "50"

extractions = {
  "Full solve [s]" : "full solve: ",
  "Iterations" : " IT=",
  "Basis setup [s]" : "Basis setup: ",
  "Krylov solve [s]" : "pCG solve: ",
  "Orthonormalization [s]" : "Gram-Schmidt: ",
  "RHS setup [s]" : "XA0X: ",
  "Source inverse [s]" : "source_inverse: ",
  "Approximation error" : "AE: ",
  "Maxnorm" : "maxnorm_final: ",
  "Onenorm" : "onenorm_final: ",
}

bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 10 ./testgeneo"

# Coarse only tests

nev =  list(range(1,26))

nev_approx_skyscrapers =  list(range(3,51))
nev_approx_layers =  list(range(2,31))
nev_approx =  nev_approx_skyscrapers
#nev = ["10", "15", "20", "25", "30"]


if True:
  permutations = OrderedDict([
    ("nev", nev_approx),#, "20", "25", "30"]),
    ("method", ["geneo", "fastrndgeneo"])
  ])

  #run(
  #  bincall = bincall,
  #  plotprefix = "Coarse approximation error with low-acc arpack",
  #  configuration = config_approx, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")

  #config_approx["coarse_reduction_test"] = "true"
  #config_approx["nev"] = "30"
  permutations = OrderedDict([
    #("nev", nev_approx),
    ("nev", list(range(1,21,3))),

    ("method", ["geneo", "geneo_1e-3", "fastrndgeneo"]), #, "fastrndgeneo2"
    ("dirichlet_interpol_in_coarse_only", ["true", "false"])
  ])

  run(
    bincall = bincall,
    plotprefix = "Coarse approximation error",
    configuration = config_approx, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")
  exit(0) #FIXME

  config_approx["overlap"] = "15"
  permutations = OrderedDict([
    ("nev", nev_approx),#, "20", "25", "30"]),
    ("method", ["geneo", "geneo_1e-3", "fastrndgeneo"]) #, "fastrndgeneo2"
  ])

  #run(
  #  bincall = bincall,
  #  plotprefix = "Coarse approximation error ovlp 15",
  #  configuration = config_approx, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")

  config_approx["part_unity"] = "sarkis"
  permutations = OrderedDict([
    ("nev", nev_approx),#, "20", "25", "30"]),
    ("method", ["geneo", "fastrndgeneo"]) #1e-3 fails! , "fastrndgeneo2"
  ])

  #run(
  #  bincall = bincall,
  #  plotprefix = "Coarse approximation error ovlp 15 Sarkis",
  #  configuration = config_approx, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")


  #exit(0)

  # Scaling per Cells

  permutations = OrderedDict([
    ("cells", ["200", "400", "600", "800"]),
    ("method", ["geneo", "geneo_1e-6", "geneo_1e-3", "fastrndgeneo", "fastrndgeneo2", "onelevel"]),
  ])

  #run(
  #  bincall = bincall,
  #  plotprefix = "Cells ",
  #  configuration = config, permutations = permutations, extractions = extractions)



  # Robustness wrt Contrast

  permutations = OrderedDict([
    ("contrast", ["1", "1e2", "1e4", "1e6"]),
    ("nev", ["10"]),#, "20", "25", "30"]),
    ("method", ["geneo", "geneo_1e-6", "geneo_1e-3", "fastrndgeneo", "fastrndgeneo2", "onelevel"]),
  ])

  #run(
  #  bincall = bincall,
  #  plotprefix = "Contrast 5x5",
  #  configuration = config, permutations = permutations, extractions = extractions)

#nev = ["5", "10", "15", "20", "30", "40", "50", "60"]
nev =  list(range(2,36,1))
#nev = ["10", "15", "20", "25", "30"]
#nev =  ["8", "10", "15", "20", "25", "30"]


# Abort criterion

config["method"] = "fastrndgeneoadaptive"
config["hybrid"] = "true"
permutations = OrderedDict([
  ("nev", nev),#, "20", "25", "30"]),
])
run(
  bincall = bincall,
  plotprefix = "Abort criterion ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")
config["hybrid"] = "false"

#exit(0) # FIXME

config["hybrid"] = "true"
permutations = OrderedDict([
  ("nev", nev),#, "20", "25", "30"]),
  ("method", ["geneo", "fastrndgeneo", "fastrndgeneoadaptive"]),
])
run(
  bincall = bincall,
  plotprefix = "Cost of abort criterion ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")
config["hybrid"] = "false"





# Scaling wrt number of cells
config["nev"] = "20"
config["hybrid"] = "true"
permutations = OrderedDict([
  ("cells", ["400", "800", "1600", "3200"]),
  ("method", ["geneo", "geneo_1e-3", "fastrndgeneo"]),
#  ("hybrid", ["true", "false"])
])
run(
  bincall = bincall,
  plotprefix = "Hybrid scaling cells ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "1/h", ylogscale = True)
config["hybrid"] = "false"



# Effectiveness/Cost of #EV

permutations = OrderedDict([
  ("nev", nev),#, "20", "25", "30"]),
  ("method", ["geneo", "geneo_1e-3", "fastrndgeneo"]),
  ("hybrid", ["true", "false"])
])
run(
  bincall = bincall,
  plotprefix = "Hybrid vs Additive EV 250 Cells ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")

permutations = OrderedDict([
  ("nev", nev),#, "20", "25", "30"]),
  ("method", ["geneo", "geneo_1e-6", "geneo_1e-3", "fastrndgeneo", "fastrndgeneo2"])
])
run(
  bincall = bincall,
  plotprefix = "EV 250 Cells ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")


config["hybrid"] = "true"
permutations = OrderedDict([
  ("nev", nev),#, "20", "25", "30"]),
  ("method", ["geneo", "geneo_1e-6", "geneo_1e-3", "fastrndgeneo", "fastrndgeneo2"])
])
run(
  bincall = bincall,
  plotprefix = "Hybrid EV 250 Cells ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")
config["hybrid"] = "false"



permutations = OrderedDict([
  ("nev", nev),#, "20", "25", "30"]),
  ("method", ["geneo", "geneo_1e-3", "fastrndgeneo", "fastrndgeneo2"]),
  ("hybrid", ["true", "false"])
])
run(
  bincall = bincall,
  plotprefix = "With fast2 Hybrid vs Additive EV 250 Cells ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")


permutations = OrderedDict([
  ("nev", list(range(10,36,1))),#, "20", "25", "30"]),
  ("method", ["geneo", "geneo_1e-3", "fastrndgeneo"]),
  ("hybrid", ["true", "false"])
])

run(
  bincall = bincall,
  plotprefix = "Many EV Hybrid vs Additive EV 250 Cells ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")


permutations = OrderedDict([
  ("nev", list(range(10,36,1))),#, "20", "25", "30"]),
  ("method", ["geneo", "geneo_1e-3", "fastrndgeneo", "fastrndgeneo2"]),
  ("hybrid", ["true", "false"])
])

run(
  bincall = bincall,
  plotprefix = "With fast2 Many EV Hybrid vs Additive EV 250 Cells ",
  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")




config["hybrid"] = "true"
permutations = OrderedDict([
  ("nev", nev),#, "20", "25", "30"]),
  ("method", ["geneo", "geneo_1e-6", "geneo_1e-3", "fastrndgeneo", "fastrndgeneo2"])
])

#run(
#  bincall = bincall,
#  plotprefix = "Hybrid EV 250 Cells ",
#  configuration = config, permutations = permutations, extractions = extractions, xlabel = "Number of eigenvectors")
config["hybrid"] = "false"
