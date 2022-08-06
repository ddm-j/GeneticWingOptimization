import numpy as np
import random
from random import randint
import time
import subprocess
import os.path
import os
import signal
import glob
from geometric_functions import *
import operator

def decision(probability):

    return random.random() < probability

#dna1 = rand_dna()
#dna2 = rand_dna()

#child1,child2 = single_point_crossover(dna1,dna2)

#print(dna1)
#print(dna2)
#print(child1)
#print(child2)

#for i in range(0,10):

#    print(decision(0.1))


class chromosome:

    def __init__(self):

        self.genes = ['sect','sect','sect','length','chord','percent_l','sweep','sweep',
                      'dihedral','dihedral','taper','taper','AoA','AoA','AoA']

        self.camber = np.arange(1, 10, 1)  # nine possibilities
        self.camber_loc = np.arange(1, 10, 1)  # nine possibilities
        self.thickness = np.arange(10, 41, 1)  # 30 possibilities
        self.percent_l = np.linspace(0.4, 0.7, 20)  # 4 possibilities
        self.sweep = np.linspace(0, 20, 40)  # 11
        self.dihedral = np.linspace(-4, 4, 20)  # 9 possibilities
        self.taper = np.linspace(0.4, 0.9, 20)  # 6 possibilities
        self.AoA = np.linspace(-4, 15, 40)  # 38 possibilities

    def rand_dna(self):

        dna_ind = str()

        for i in range(0,3):

            dna_ind += str(random.choice(range(len(self.camber)))) + str(random.choice(range(len(self.camber_loc))))\
                       + str(random.choice(range(len(self.thickness)))) + '-'

        sect1 = str(random.choice(self.camber)) + str(random.choice(self.camber_loc)) + str(random.choice(self.thickness))
        sect2 = str(random.choice(self.camber)) + str(random.choice(self.camber_loc)) + str(random.choice(self.thickness))
        sect3 = str(random.choice(self.camber)) + str(random.choice(self.camber_loc)) + str(random.choice(self.thickness))

        taper1 = random.choice(range(len(self.taper)))

        dna_ind += str(random.choice(range(len(self.percent_l)))) + '-'
        dna_ind += str(random.choice(range(len(self.sweep)))) + '-'
        dna_ind += str(random.choice(range(len(self.sweep)))) + '-'
        dna_ind += str(random.choice(range(len(self.dihedral)))) + '-'
        dna_ind += str(random.choice(range(len(self.dihedral)))) + '-'

        taper_list = [x for x in self.taper if x < self.taper[taper1]]

        if taper_list == []:

            taper2 = taper1

        else:

            taper2 = random.choice(range(len(taper_list)))

        dna_ind += str(taper1) + '-'
        dna_ind += str(taper2) + '-'
        dna_ind += str(random.choice(range(len(self.AoA)))) + '-'
        dna_ind += str(random.choice(range(len(self.AoA)))) + '-'
        dna_ind += str(random.choice(range(len(self.AoA))))

        self.dna = [sect1, sect2, sect3, 20, 4, random.choice(self.percent_l), random.choice(self.sweep),
               random.choice(self.sweep), random.choice(self.dihedral), random.choice(self.dihedral),
               taper1, taper2, random.choice(self.AoA), random.choice(self.AoA), random.choice(self.AoA)]

        self.dna_ind = dna_ind

    def mutate_gene(self,gene_id):

        if gene_id == 'sect':

            mutated = str(random.choice(self.camber)) + str(random.choice(self.camber_loc)) + str(random.choice(self.thickness))

            return mutated

        elif gene_id == 'length':

            return 20

        elif gene_id == 'chord':

            return 4

        else:

            mutated = random.choice(getattr(self,gene_id))

            return mutated




class individual:

    def __init__(self,gen_id,case,dna):

        self.sys_path = '/home/bjorn/OpenFOAM/bjorn-5.0/run/genetic_cases/gen_'+str(gen_id)+'/'
        self.case = case
        self.folder = '/postProcessing/forces/0/forces.dat'
        self.dna = dna
        self.file_path = self.sys_path + self.case + self.folder
        self.lift = np.NaN
        self.drag = np.NaN
        self.moment = np.NaN

    def get_random_dna(self):

        chrom = chromosome()
        chrom.rand_dna()
        self.dna = chrom.dna


    def line2dict(self,line):
        tokens_unprocessed = line.split()
        tokens = [x.replace(")", "").replace("(", "") for x in tokens_unprocessed]
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

    def follow(self,file):

        file.seek(0, 2)
        while True:
            line = file.readline()
            if not line:
                time.sleep(0.1)
                continue
            yield line

    def monitor_sim(self):

        datafile = open(self.file_path, 'r')
        datalines = self.follow(datafile)

        t = []
        drag = []
        lift = []
        moment = []

        i = 0

        conv = 100

        for line in datalines:
            print("New Line")
            data_dict = self.line2dict(line)
            t += [data_dict['time']]
            drag += [data_dict['force']['pressure'][0] + data_dict['force']['viscous'][0]]
            lift += [data_dict['force']['pressure'][1] + data_dict['force']['viscous'][1]]
            moment += [data_dict['moment']['pressure'][2] + data_dict['moment']['viscous'][2]]

            if i == 0:
                i += 1
                continue

            aero = np.array(lift) / np.array(drag)

            if i > 10:
                conv = 100 * abs((aero[-1] - aero[-10]) / aero[-10])

                print('Current convergence:'+str(conv)+'%')

            if conv < 5.0 and i > 10:
                return np.mean(lift[-5:]), np.mean(drag[-5:]), np.mean(moment[-5:])

            elif i > 45:

                return np.mean(lift[-5:]), np.mean(drag[-5:]), np.mean(moment[-5:])

            i += 1

    def run_sim(self,genMesh=True):

        os.chdir(self.sys_path+self.case)

        #Generate Mesh

        if genMesh == True:

            p = subprocess.Popen(['/bin/sh','-c','./genMesh'],cwd=self.sys_path+self.case)#,stdout=FNULL, stderr=subprocess.STDOUT)
            p.wait()

        #Run Simulation

        p = subprocess.Popen(['/bin/sh','-c','./runSim'],cwd=self.sys_path+self.case, preexec_fn=os.setsid)

        while not(os.path.isfile(self.file_path)):

            time.sleep(1)

        self.lift,self.drag,self.moment = self.monitor_sim()

        os.killpg(os.getpgid(p.pid),signal.SIGTERM)

        # Reconstruct Volumetric Fields

        p = subprocess.Popen(['/bin/sh','-c','./reconstructFields'],cwd=self.sys_path+self.case)

        p.wait()



class generation:

    os.chdir('/home/bjorn/OpenFOAM/bjorn-5.0/run/genetic_cases/')

    def __init__(self,id):

        self.gen_id = id
        self.dir = '/home/bjorn/OpenFOAM/bjorn-5.0/run/genetic_cases/'
        self.gen_dir = self.dir + 'gen_' + str(self.gen_id) + '/'
        self.individuals = []
        self.results = {}

    def rand_gen(self,n_individuals):

        gen_list = []

        for i in range(0,n_individuals):

            tmp = individual(self.gen_id,'indiv_'+str(i+1),0)
            tmp.get_random_dna()

            gen_list = np.append(gen_list,tmp)

        self.individuals = gen_list

    def read_gen(self,children):

        gen_list = []

        for i in range(0, len(children)):
            tmp = individual(self.gen_id, 'indiv_' + str(i + 1), children[i])

            gen_list = np.append(gen_list, tmp)

        self.individuals = gen_list

    def single_point_crossover(self, dna1, dna2):

        point = randint(1, 15)

        child1 = dna1[:point] + dna2[point:]
        child2 = dna2[:point] + dna1[point:]

        return child1, child2

    def case_files(self):

        p = subprocess.Popen(['mkdir','-p','gen_'+str(self.gen_id)],cwd=self.dir)

        for i in range(0,len(self.individuals)):

            p = subprocess.Popen(['mkdir','-p','indiv_'+str(i+1)],cwd=self.gen_dir)
            copy_folders = glob.glob(self.dir+'base_case/*')
            for j in range(0,len(copy_folders)):
                p = subprocess.Popen(['cp','-rf',copy_folders[j],self.gen_dir+self.individuals[i].case+'/'],cwd=self.dir)

            wing_design_pipeline(self.individuals[i].dna,[0.3,0.7,1.4,2.5],self.gen_dir+self.individuals[i].case+'/constant/triSurface/')


    def gen_sim(self):

        for i in range(0,len(self.individuals)):

            self.individuals[i].run_sim()
            self.results[self.individuals[i].case] = [self.individuals[i].dna,self.individuals[i].lift/self.individuals[i].drag]

    def selection(self,n_parents):

        total_fitness = 0
        p = {}

        for indiv_id in self.results:

            total_fitness += self.results[indiv_id][1]

        p_prev = 0

        for indiv_id in self.results:

            p_curr = self.results[indiv_id][1]/total_fitness + p_prev

            p[indiv_id] = round(p_curr,4)

            p_prev = p_curr

        sorted_p = sorted(p.items(), key=operator.itemgetter(1))

        parents = []

        for i in range(0,n_parents):

            probability = random.uniform(0,1)

            for j in range(0,len(sorted_p)):

                if probability < sorted_p[j][1]:

                    parents += [sorted_p[j][0]]

                    break

        return parents

    def mate_parents(self,parents):

        children = []

        for i in range(0,len(parents),2):

            parent1 = self.results[parents[i]][0]
            parent2 = self.results[parents[i+1]][0]

            child1,child2 = self.single_point_crossover(parent1,parent2)

            children += [child1]
            children += [child2]

        return children


    def mutate(self,children,rate=0.01):

        tmp_chromo = chromosome()

        mutated_children = children.copy()

        for i in range(0,len(children)):

            for j in range(0,len(children[0])):

                p = random.uniform(0,1)

                if p < rate:

                    # Mutate

                    mutated = tmp_chromo.mutate_gene(tmp_chromo.genes[j])

                    mutated_children[i][j] = mutated

        return mutated_children


    def test_run(self):

        for i in range(0, len(self.individuals)):

            lift = (i+1)*10.0
            drag = (i+2)*5.0

            self.results[self.individuals[i].case] = [self.individuals[i].dna,
                                                  lift / drag]


    def test_write(self):

        fh = open('genResults.txt','w')

        lines = []

        for i in range(0,len(self.results)):

            lines = np.append(lines,self.results[i])

        fh.writelines(lines)

        fh.close()


class evolutionary_simulation:

    def __init__(self,n_generations,n_population,mut_rate):

        self.n_gen = n_generations
        self.n_pop = n_population
        self.mutation_rate = mut_rate


    def run(self):

        self.gens = []

        for i in range(0,self.n_gen):

            # Create Generation

            print('Beginning Generation: ',i)

            tmp_gen = generation(i)

            gen_list = []

            # If i == 0, Generation is Random

            if i == 0:

                tmp_gen.rand_gen(self.n_pop)

            else:

                tmp_gen.read_gen(tmp_children)

            # Create Case Files For Generation

            tmp_gen.case_files()

            # Run Simulation for Generation

            tmp_gen.gen_sim()

            # Select Parents Based on Results

            tmp_parents = tmp_gen.selection(self.n_pop)

            # Mate Parents

            tmp_children = tmp_gen.mate_parents(tmp_parents)

            # Mutate Children

            tmp_children = tmp_gen.mutate(tmp_children,self.mutation_rate)

            # Append Generation

            self.gens += [tmp_gen]




