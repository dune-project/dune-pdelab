from pandas import *
import matplotlib.pyplot as plt
import itertools
import os
import subprocess
import re
from collections import OrderedDict
from slurmpy import Slurm


def run (bincall, plotprefix, configuration, permutations, extractions):

  results = {}
  for tpl in list(itertools.product(*permutations.values())):
    results[tpl[1:]] = {}
    for extraction in extractions.keys():
      results[tpl[1:]][extraction] = {}


  for tpl in list(itertools.product(*permutations.values())):

    print ("Running " + str(tpl))

    runconfig = configuration.copy()
    runconfig.update(zip(permutations.keys(), tpl))
    params = ' '.join(['-%s %s' % (key, value) for (key, value) in runconfig.items()])


    #sbatch --out "test.out" batch.sbatch "params..."
    os.system("sbatch --job-name=prmstudy --out \"" + str(tpl) + ".out\" batch.sbatch \"" + params + "\"")

    #os.system (bincall + " " + params + " > out")


  print("waiting for runs...")
  os.system("while squeue | grep --quiet prmstudy; do sleep 3; done")

  for tpl in list(itertools.product(*permutations.values())):
    print ("Reading " + str(tpl) + ".out")

    with open(str(tpl) + ".out", encoding="utf-8") as f:
      for line in f:
        print(line)
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
    plt.title(extraction)
    plt.legend()
    fig.savefig(plotprefix + extraction + ".svg")
    fig.clf()
  #plt.show()











config = {
  "cells": "200",
  "subdomainsx": "4",
  "subdomainsy": "1",
  "nev": "20",
  "nev_arpack": "20",
  "overlap": "1",
  "layers": "50",
  "contrast": "1e6",
  "method": "geneo",
  "part_unity": "standard",
  "hybrid": "false"
}

extractions = {
  "Full solve" : "full solve: ",
  "Iterations" : " IT=",
  "Basis setup" : "Basis setup: ",
  "pCG solve" : "pCG solve: "
}



# Scaling per Cells

permutations = OrderedDict([
  ("cells", ["200", "400", "600", "800"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo", "fastrndgeneo2", "onelevel"]),
])

#run(
#  bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 4 ./testgeneo",
#  plotprefix = "Cells ",
#  configuration = config, permutations = permutations, extractions = extractions)



# Robustness wrt Contrast

permutations = OrderedDict([
  ("contrast", ["1", "1e2", "1e4", "1e6"]),
  ("nev", ["10"]),#, "20", "25", "30"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo", "fastrndgeneo2", "onelevel"]),
])

#run(
#  bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 25 ./testgeneo",
#  plotprefix = "Contrast 5x5",
#  configuration = config, permutations = permutations, extractions = extractions)


# Effectiveness/Cost of #EV


permutations = OrderedDict([
  ("nev", ["5", "10", "15", "20", "30", "40", "50", "60"]),#, "20", "25", "30"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo", "fastrndgeneo2"]),
  ("cells", ["100"])
])

run(
  bincall = "--export=OMP_NUM_THREADS=1 ./testgeneo",
  plotprefix = "EV 100 Cells ",
  configuration = config, permutations = permutations, extractions = extractions)

exit(0)


permutations = OrderedDict([
  ("nev", ["5", "10", "15", "20", "30", "40", "50", "60"]),#, "20", "25", "30"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo", "fastrndgeneo2"]),
  ("cells", ["250"])
])

#run(
#  bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 4 ./testgeneo",
#  plotprefix = "EV 250 Cells ",
#  configuration = config, permutations = permutations, extractions = extractions)


permutations = OrderedDict([
  ("nev", ["5", "10", "15", "20", "30", "40", "50", "60"]),#, "20", "25", "30"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo", "fastrndgeneo2"]),
  ("cells", ["500"])
])

run(
  bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 4 ./testgeneo",
  plotprefix = "EV 500 Cells ",
  configuration = config, permutations = permutations, extractions = extractions)

exit(0) # FIXME

permutations = OrderedDict([
  ("nev", ["5", "10", "15", "20", "30", "40"]),#, "20", "25", "30"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo", "fastrndgeneo2"]),
  ("cells", ["1000"])
])

run(
  bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 4 ./testgeneo",
  plotprefix = "EV 1000 Cells ",
  configuration = config, permutations = permutations, extractions = extractions)













permutations = OrderedDict([
  ("nev", ["5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "20", "25", "30"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo"]),
  ("cells", ["800"]),
])

#run(
#  bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 4 ./testgeneo",
#  plotprefix = "EV 800 Cells ",
#  configuration = config,
#  permutations = permutations,
#  extractions = extractions
#  )

permutations = OrderedDict([
  ("nev", ["5", "10", "15", "20", "25", "30"]),
  ("method", ["geneo", "inexactgeneo", "fastrndgeneo"]),
  ("cells", ["800"]),
  ("subdomainsx", ["4"]),
  ("subdomainsy", ["4"]),
])

#run(
#  bincall = "OMP_NUM_THREADS=1 && export OMP_NUM_THREADS && make testgeneo && mpirun -x OMP_NUM_THREADS -np 16 ./testgeneo",
#  plotprefix = "EV 800 Cells 4x4 ",
#  configuration = config,
#  permutations = permutations,
#  extractions = extractions
#  )
