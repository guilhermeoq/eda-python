#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rosetta import *
from toolbox import cleanATOM
from math import sqrt
import sys
from random import randrange, choice
from copy import deepcopy
from random import uniform, choice, randint, gauss

conformations = []
sequence = []
list_atoms = []
energies = []
selection = []
mean_sdev_phi = []
mean_sdev_psi = []
interval = []
f = open("saida.csv", "w")

class Aminoacid:
    def __init__(self,phi,psi):
        self.phi = phi
        self.psi = psi

class Atomo:
    def __init__(self, x, y, z,tipo):
        self.x = x
        self.y = y
        self.z = z
        self.tipo = tipo

def create_sequence():
    for i in range(p.total_residue()):
        sequence.append(Aminoacid(p.phi(i+1), p.psi(i+1)))
    return

def create_conformations(sequence):
    for k in range(100):
        sec_copy = deepcopy(sequence)
        i=0
        while i < len(sequence):
            sec_copy[i].phi = uniform(-180,180) #(-30,-120)
            sec_copy[i].psi = uniform(-180,180) #(-50,0)
            i+=1
        conformations.append(sec_copy)
    return

def create_conformations_by_model(sequence):
    for i in range(len(conformations)):
        for k in range(len(conformations[i])):
            conformations[i][k].phi = gauss(mean_sdev_phi[k][0], mean_sdev_phi[k][1])
            conformations[i][k].psi = gauss(mean_sdev_psi[k][0], mean_sdev_psi[k][1])
    return

def calc_energies():
    energies[:] = []
    size_seq = len(sequence)
    size_conf = len(conformations)
    for c in range(size_conf):
        for r in range(p.total_residue()):
            p.set_phi(r+1,float( conformations[c][r].phi ))
            p.set_psi(r+1,float( conformations[c][r].psi ))
        xyz = extract_coordinates_from_pose_3x1(p)
        natoms = len(xyz[0])
        for i in range(natoms):  #in range(390):
            list_atoms[i].x = xyz[0][i]
            #print "X: "+str(list_atoms[i].x)
            list_atoms[i].y = xyz[1][i]
            #print "Y: "+str(list_atoms[i].y)
            list_atoms[i].z = xyz[2][i]
            #print "Z: "+str(list_atoms[i].z)
        top = iterate_vdw(len(list_atoms))
        #print "Apenas VDW: " + str(top)
        #top = top + iterate_ech(len(list_atoms))
        #print "VDW + ECH: " + str(top)
        new_element = (top,c)
        #print "Calculated: " + str(new_element[0])
        energies.append(new_element)
    analysis()
    return

def tournment():
    size_conf = len(conformations)
    num_selections = int((size_conf * 400)/100)
    i=0
    size_conf -= 1
    while i < num_selections:
        player1 = randint(1,size_conf)
        player2 = randint(1,size_conf)
        winner,loser = play_duel(player1,player2)
        selection.append(energies[winner])
        i+=1
    return

def cal_model():
    size_sequence = len(sequence)
    size_selection = len(selection)
    #size_conformations = len(conformations)

    for k in range(size_sequence):
        sum_phi=0
        sum_psi=0
        
        for i in range(size_selection):
            sum_phi += conformations[selection[i][1]][k].phi
            sum_psi += conformations[selection[i][1]][k].psi

        mean_phi = sum_phi / size_selection
        mean_psi = sum_psi / size_selection

        variance_phi = 0
        variance_psi = 0

        for i in range(size_selection):
            variance_phi += ((conformations[selection[i][1]][k].phi - mean_phi)**2)/size_selection
            variance_psi += ((conformations[selection[i][1]][k].psi - mean_psi)**2)/size_selection

        standard_deviation_phi = sqrt(variance_phi)
        standard_deviation_psi = sqrt(variance_psi)

        value_phi = (mean_phi, standard_deviation_phi)
        value_psi = (mean_psi, standard_deviation_psi)
        
        mean_sdev_phi.append(value_phi)
        mean_sdev_psi.append(value_psi)
    return

def play_duel(player1,player2):
    if (energies[player1][0] <= energies[player2][0]):
        winner, loser = player1, player2
    else:
        winner, loser = player2, player1
    return winner, loser

def create_interval():
    for i in range(41):
        interval.append(i*500)
    return

def calc_mean_desv_population():
    size_energies = len(energies)
    #size_conformations = len(conformations)
    sum_energies = 0
    for i in range(size_energies):
        sum_energies += energies[i][0]

    mean_energies = sum_energies / size_energies

    variance_energies = 0

    for i in range(size_energies):
        variance_energies += ((energies[i][0] - mean_energies)**2)/size_energies
    
    standard_deviation_energies = sqrt(variance_energies)
    value_energies = (mean_energies, standard_deviation_energies)
    
    f.write("Mean and Std. Dev. Energies of Population: "+str(value_energies)+"\n")
    return

def analysis():
    i = 1
    while i < len(interval):
        count=0
        for j in range(len(energies)):
            if (energies[j][0] < interval[i]) and (energies[j][0] >= interval[i-1]):
                count+=1
        print("Interval: "+str(interval[i-1])+" - "+str(interval[i])+ "; Count: "+str(count))
        f.write(str(interval[i-1])+" - "+str(interval[i])+";"+str(count)+"\n")
        i+=1
    f.write("end\n")
    calc_mean_desv_population()
    return

def execute():
    for i in range(5):
        print "> Calculating energies...\n"
        calc_energies()
        print "> Realizing a tournment...\n"
        tournment()
        print "> Calculating Mean and Standard Deviation...\n"
        cal_model()
        print "> Generating new population by model...\n"
        create_conformations_by_model(sequence)
        #print_file(i)
    return

def open_file(nome_arquivo):
    i=0
    arquivo = open(str(nome_arquivo),"r")
    linha = arquivo.readline()
    header = linha.split()
    num_atomos = int(header[0])

    while i < num_atomos:
        linha = arquivo.readline()
        linha = linha.split()
        tipo=linha[1]
        list_atoms.append(Atomo(float(linha[2]),float(linha[3]),float(linha[4]),tipo[0]))
        i+=1
    arquivo.close()
    return

'''Obter raio de van der Waals'''
def get_radii(atom):
    if atom.tipo == 'C':
        return 1.7
    elif atom.tipo == 'H':
        return  1.2
    elif atom.tipo == 'N':
        return 1.55
    elif atom.tipo == 'O':
        return 1.52
    elif atom.tipo == 'S':
        return 1.8
    else:
        return 0.0

def get_charge(atom):
    if atom.tipo == 'C':
        return 12.011
    elif atom.tipo == 'H':
        return  1.008
    elif atom.tipo == 'N':
        return 14.007
    elif atom.tipo == 'O':
        return 15.999
    elif atom.tipo == 'S':
        return 32.060
    else:
        return 0.0

'''Calcula a Distancia Euclidiana entre dois atomos'''
def dist_euclid(atom1,atom2):
    result = sqrt(((atom1.x - atom2.x) * (atom1.x - atom2.x)) +\
        ((atom1.y - atom2.y) * (atom1.y - atom2.y)) +\
        ((atom1.z - atom2.z) * (atom1.z - atom2.z)))
    return result

'''Calcula van der Waals entre dois atomos '''
def calc_vdw(atom1,atom2):
    d = dist_euclid(atom1,atom2)
    r = d / ( get_radii(atom1) + get_radii(atom2) )

    r6 = r*r*r*r*r*r
    r12 = r6*r6
    result = 0

    if r <= 0.8: #if r <= 0.80:
        r=0.8 #r=0.80
        r6  = r*r*r*r*r*r
        r12 = r6*r6

        result = (1/r12) - (2/r6)

    elif (r > 0.8) and (d < 8.0): #elif (r > 0.8) and (d < 8.0):
        r6  = r*r*r*r*r*r
        r12 = r6*r6
        result = (1/r12) - (2/r6)
    return result

def calc_ech(atom1,atom2):
    dieletric = 1
    d = dist_euclid(atom1,atom2)
    r = d / ( get_radii(atom1) + get_radii(atom2) )
    result = 0.0

    if r <= 0.75:
        r = 0.75
        result = (get_charge(atom1) * get_charge(atom2))/(dieletric * r)

    elif (r > 0.75) and (d < 13.0):
        result = (get_charge(atom1) * get_charge(atom2))/(dieletric * r)
    return result

def iterate_vdw(num_atomos):
    k = num_atomos-1
    n = num_atomos
    potencial=0.0

    for i in range(n-2):
        j=i+1
        while j < n:
            potencial+=calc_vdw(list_atoms[i],list_atoms[j])
            j += 1
    return potencial

def iterate_ech(num_atomos):
    k = num_atomos-1
    n = num_atomos
    potencial=0.0

    for i in range(n-2):
        j=i+1
        while j < n:
            potencial+=calc_ech(list_atoms[i],list_atoms[j])
            j += 1
    return potencial

def get_angles():
	for i in range(p.total_residue()):
		phis = p.phi(i + 1)
		psis = p.psi(i + 1)
		#chis = p.residue(i+1).chi()


def random_structure():
	for i in range(p.total_residue()):
		p.set_phi(i+1,float(choice(range(-180,180,1))))
		p.set_psi(i+1,float(choice(range(-180,180,1))))

def print_file(num):
	p.dump_pdb("output_"+str(num)+".pdb")

def extract_coordinates_from_pose_3x1( pose , selection = [] ,
        atom_names = [] , not_atom_names = [] ):
    # default to all residues
    if not selection:
        selection = range( 1 , pose.total_residue() + 1 )
    # empty coordinate holders
    coords = [[] for i in range(3)]
    # for each residue
    for resi in selection:
        resi = pose.residue(resi)
        # for each atom
        for atom in range( resi.natoms() ):
            # check if the atom is a type desired
            if ( not atom_names or resi.atom_name( atom + 1 ).strip()
                    in atom_names ) and ( not
                    resi.atom_name( atom + 1 ).strip() in not_atom_names ):
                atom = resi.xyz( atom + 1 )
                # store the coordinates
                coords[0].append( atom[0] )    # x
                coords[1].append( atom[1] )    # y
                coords[2].append( atom[2] )    # z
    return coords

'''Main'''
if __name__ == "__main__":
    print "Iniciando Rosetta...\n"
    rosetta.init()
    p = Pose()
    print "> Opening the PDB file...\n"
    pose_from_pdb(p,"1A11.clean.pdb")
    print "> Opening the XYZ file...\n"
    open_file(str("1A11.xyz"))
    print "> Identifying dihedral angles...\n"
    create_sequence()
    print "> Creating random conformations...\n"
    create_conformations(sequence)
    print "> Creating interval for analysis...\n"
    create_interval()
    execute()