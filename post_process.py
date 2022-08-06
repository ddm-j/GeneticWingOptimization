import numpy as np
import matplotlib.pyplot as plt

def line2dict(line):
    tokens_unprocessed = line.split()
    tokens = [x.replace(")","").replace("(","") for x in tokens_unprocessed]
    floats = [float(x) for x in tokens]
    data_dict = {}
    data_dict['time'] = floats[0]
    force_dict = {}
    force_dict['pressure'] = floats[1:4]
    force_dict['viscous'] = floats[4:7]
    force_dict['porous'] = floats[7:10]
    moment_dict = {}
    moment_dict['pressure'] = floats[10:13]
    moment_dict['viscous'] = floats[13:16]
    moment_dict['porous'] = floats[16:19]
    data_dict['force'] = force_dict
    data_dict['moment'] = moment_dict
    return data_dict


forces_file = '/home/bjorn/OpenFOAM/bjorn-5.0/run/genetic_cases/gen_1/indiv_4/postProcessing/forces/0/forces.dat'

datafile = open(forces_file, 'r')

t = []
drag = []
lift = []
moment = []

for line in datafile:

    if line[0] == '#':
        continue

    data_dict = line2dict(line)
    t += [data_dict['time']]
    drag += [data_dict['force']['pressure'][0] + data_dict['force']['viscous'][0]]
    lift += [data_dict['force']['pressure'][1] + data_dict['force']['viscous'][1]]
    moment += [data_dict['moment']['pressure'][2] + data_dict['moment']['viscous'][2]]

print(lift)

aero = np.array(lift)/np.array(drag)

plt.plot(aero)
plt.title('Wing Proto - Aerodynamic Coefficient Convergence')
plt.xlabel('Iterations (x10)')
plt.ylabel('Lift Over Drag')
plt.show()
