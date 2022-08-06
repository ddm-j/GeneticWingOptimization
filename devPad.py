import numpy as np
from time import time
import random
from geometric_functions import *
from genetic_functions import *
import matplotlib.pyplot as plt

forces_file = '/home/bjorn/OpenFOAM/bjorn-5.0/run/genetic_steady/postProcessing/forces/0/forces.dat'

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


def follow(file):

    file.seek(0,2)
    while True:
        line = file.readline()
        if not line:
            time.sleep(0.1)
            continue
        yield line


def monitor_sim():

    datafile = open(forces_file,'r')
    datalines = follow(datafile)

    t = []
    drag = []
    lift = []
    moment = []

    i = 0

    for line in datalines:
        print('New Line')
        data_dict = line2dict(line)
        t += [data_dict['time']]
        drag += [data_dict['force']['pressure'][0] + data_dict['force']['viscous'][0]]
        lift += [data_dict['force']['pressure'][1] + data_dict['force']['viscous'][1]]
        moment += [data_dict['moment']['pressure'][2] + data_dict['moment']['viscous'][2]]

        if i == 0:
            i+=1
            continue

        aero = np.array(lift)/np.array(drag)

        conv = 1000

        if i > 11:

            conv = 100*abs(aero[-1]-aero[-10])/aero[-10]

            print(conv)

        if conv < 1.0 and i > 10:

            return lift[-5:].mean(),drag[-5:].mean(),moment[-5:].mean()

        i += 1

lift,drag,moment = monitor_sim()

print('Simulation Convergence Reached.')
print('Lift',lift,'Drag',drag,'L/D',lift/drag)