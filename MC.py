'''
Descripttion: 
version: 
Author: kalice
Date: 2021-03-08 08:48:05
LastEditors: kalice
LastEditTime: 2021-03-18 23:15:05
'''
# coding=utf-8
import ase
import ase.io
import glob
import random
import numpy as np
import time
import math

lmp = 'conf.lmp'
atoms_total = ase.io.read(lmp, format='lammps-data', Z_of_type=[26, 15, 8, 3], style="atomic")#读取lammps模拟数据


def get_elements_position(atoms, Z_of_type):#读取原子结构参数
    list_elements_position = []
    for i, num in enumerate(atoms.numbers):
        if num == Z_of_type:
            list_elements_position.append(list(atoms.positions[i]))
    return list_elements_position


def get_elements_atoms(atoms, Z_of_type):
    elements_indx = []
    for i, num in enumerate(atoms.numbers):
        if num == Z_of_type:
            elements_indx.append(i)
    atoms_elements = atoms[elements_indx]
    return atoms_elements


def delete_Li_atom(atoms, num):
    Li_indx = []
    atoms_copy = atoms.copy()
    for i, Z in enumerate(atoms.numbers):
        if Z == 3:
            Li_indx.append(i)
    del atoms_copy[Li_indx]
    random.shuffle(Li_indx)
    for ii in range(num):
        atoms_copy.append(atoms[Li_indx[ii]])
    return atoms_copy


def MC(energy1, energy2, pop1, pop2):#简单计算锂离子玻尔兹曼分布
    p = min(1, math.exp(-(energy2-energy1)/0.0256))
    if p >= random.random():
        return energy2, pop2
    else:
        return energy1, pop1



from deepmd.calculator import DP #调用deepmd进行能量计算
random.seed(time.time())

#对lammps数据进行预处理
Li_position = get_elements_position(atoms_total, 3)
atoms_Li = get_elements_atoms(atoms_total, 3)
atoms_framework = delete_Li_atom(atoms_total, 0)

atoms_target = delete_Li_atom(atoms_total, 50)
Li_position_target = get_elements_position(atoms_target, 3)
# Get init energy #
atoms_target.calc = DP(model="frozen_model.pb")
init_energy = atoms_target.get_potential_energy()[0]
# Get DNA of target #
pop = np.array([1 if ii in Li_position_target else 0 for ii in Li_position])
# Set numbers of MC #
num_MC = 10
# MC process #
txt = open('result.txt', 'w+')
txt.write(str(time.time()) + '\n')
txt.write('Num  Energy/eV  Pop\n')
for i in range(num_MC):
    txt.write('{}  {}  {}\n'.format(i, init_energy, pop))
    pop_copy = pop.copy()
    atoms_framework_copy = atoms_framework.copy()
    random.shuffle(pop_copy)
    atoms_framework_copy.extend(atoms_Li[pop_copy.astype(np.bool)])
    atoms_new = atoms_framework_copy
    # print(atoms_new.symbols)
    atoms_new.set_calculator(DP(model="frozen_model.pb"))
    energy_next = atoms_new.get_potential_energy()[0]
    init_energy, pop = MC(init_energy, energy_next, pop, pop_copy)


txt.write(str(time.time()) + '\n')
txt.close()


