from pandas import *
import matplotlib.pyplot as plt
import itertools
import os
import subprocess
import re
from collections import OrderedDict
import csv

def run (bincall, plotprefix, configuration, permutations, extractions):

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

    #os.system("sbatch --job-name=prmstudy --out \"out/" + str(tpl) + ".out\" batch.sbatch \"" + params + "\"")
    print(bincall + " " + params + " > out/" + str(tpl) + ".out")
    os.system(bincall + " " + params + " > out/\"" + str(tpl) + ".out\"")



  print("waiting for runs...")
  #os.system("while squeue | grep --quiet prmstudy; do sleep 3; done")

  for tpl in list(itertools.product(*permutations.values())):
    print ("Reading out/" + str(tpl) + ".out")

    with open("out/" + str(tpl) + ".out", encoding="utf-8") as f:
      for line in f:
        #print(line)
        for extraction in extractions.keys():
          match = re.search('(?<=' + extractions[extraction] + ').*', line)
          if (match):
            results[tpl[1:]][extraction][tpl[0]] = match.group(0)



  primaries = list(permutations.values())[0]

  for extraction in extractions.keys():

    fig = plt.figure()

    print (permutations.values())
    secondaries = list(permutations.values())
    secondaries.pop(0)
    for subtpl in list(itertools.product(*secondaries)):
      line, = plt.plot(np.array(primaries), [float(results[subtpl][extraction].get(x, 0.0)) for x in primaries])
      line.set_label(" ".join(str(subtpl)))

    plt.title(extraction)
    plt.legend()
    fig.savefig(plotprefix + extraction + ".svg")
    fig.clf()

    with open(plotprefix + extraction + ".csv", mode='w') as csv_file:
      csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
      csv_writer.writerow([list(permutations.keys())[0]] + primaries)
      for subtpl in list(itertools.product(*secondaries)):
        csv_writer.writerow([str(subtpl)] + [float(results[subtpl][extraction].get(x, 0.0)) for x in primaries])









bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 25 ./testgeneo"

extractions = {
#  "Condition" : "Condition estimate: ",
  "Iterations" : " IT=",
  "Matrix setup" : " M=",
  "Basis setup" : " G=",
  "CG Solver" : " S=",
  "Solver total" : " F="
}

config = {
  "layers": "25",
  "cells": 125,

  "overlap": "1",
  "nev": "5",
  "nev_arpack": "10",
  "use_threshold": "false",

  "layer_model": "true",
  "contrast": "1e6",

  "extend_domain_with_procs": "false",

#  "condition_estimate": "true"
  "condition_estimate": "false"
}



# Condition vs Contrast

# Layers
permutations = OrderedDict([
  ("contrast", ["1e1", "1e2", "1e3", "1e4", "1e5", "1e6"]),
  ("nev", ["-1"] + list(range(3,10)))
])
run(
  bincall = bincall,
  plotprefix = "Layer model ",
  configuration = config, permutations = permutations, extractions = extractions)



# Skyscrapers
config["layer_model"] = "false"


nev =  list(range(3,10))
permutations = OrderedDict([
  ("contrast", ["1e1", "1e2", "1e3", "1e4", "1e5", "1e6"]),
  ("nev", ["-1"] + nev)
])
run(
  bincall = bincall,
  plotprefix = "Skyscraper model ",
  configuration = config, permutations = permutations, extractions = extractions)
