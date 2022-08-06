import numpy as np
from geometric_functions import *
from genetic_functions import *

simulation = evolutionary_simulation(2,4,0.01)

simulation.run()

for i in range(0,len(simulation.gens)):

    print(simulation.gens[i].results)