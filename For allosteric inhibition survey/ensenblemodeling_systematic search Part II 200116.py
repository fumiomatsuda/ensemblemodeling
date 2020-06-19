#-------------------------------------------------------------------------------
# Name:        EnsembleModeling.py
# Purpose:
#
# Author:      Fumio_Matsuda, Osaka Univeristy
# Version: 0.10
# Created:     25/5/2018
# Copyright:   (c) Fumio_Matsuda 2016
# Licence:     CC BY SA
#-------------------------------------------------------------------------------
import re, csv, numpy, copy, time
import scipy.integrate
from math import log, log10, exp, pow
import pp
import csv



#
# Rate equations
#
def func_constant(conc, parameters):
    # for constant reaction
    Vmax = parameters["Vmax"]
    S = conc[0]
    v = Vmax
    return v
def func_irreversible_Michaelis_Menten(conc, parameters):
    # Substrate
    # Product 1
    # irreversible
    Km = parameters["Km"]
    Vmax = parameters["Vmax"]
    S = conc[0]
    v = Vmax * S /(Km + S)
    return v
def func_irreversible_Michaelis_Menten_bi_bi(conc, parameters):
    # Substrate 2
    # Product 2
    # reversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    v = Vmax * (a*b/(Kma*Kmb))/ (1.0 + a/Kma+ b/Kmb + a/Kma * b/Kmb)
    return v
def func_irreversible_Michaelis_Menten_bi_uni(conc, parameters):
    # Substrate 2
    # Product 1
    # irreversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Vmax = parameters["Vmax"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    v = Vmax * (a*b/(Kma*Kmb))/ (1.0 + a/Kma + b/Kmb + a/Kma * b/Kmb)
    return v

def func_irreversible_Michaelis_Menten_i(conc, parameters):
    # Substrate
    # Product 1
    # irreversible
    Km = parameters["Km"]
    Vmax = parameters["Vmax"]
    S = conc[0]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))* Vmax * S /(Km + S)
    return v
def func_irreversible_Michaelis_Menten_bi_bi_i(conc, parameters):
    # Substrate 2
    # Product 2
    # reversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]

    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))* Vmax * (a*b/(Kma*Kmb))/ (1.0 + a/Kma+ b/Kmb + a/Kma * b/Kmb)
    return v
def func_irreversible_Michaelis_Menten_bi_uni_i(conc, parameters):
    # Substrate 2
    # Product 1
    # irreversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Vmax = parameters["Vmax"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))* Vmax * (a*b/(Kma*Kmb))/ (1.0 + a/Kma + b/Kmb + a/Kma * b/Kmb)
    return v

def func_mass_action_three_substrate(conc, parameters):
    # Substrate 3
    # Product nd
    # irreversible
    Vmax = parameters["Vmax"]
    S0 = conc[0]
    S1 = conc[1]
    S2 = conc[2]
    v = Vmax * S0 * S1 * S2
    return v
def func_mass_action_two_substrate(conc, parameters):
    # Substrate 2
    # Product nd
    # irreversible
    Vmax = parameters["Vmax"]
    S0 = conc[0]
    S1 = conc[1]
    v = Vmax * S0 * S1
    return v
def func_mass_action_one_substrate(conc, parameters):
    # Substrate 1
    # Product nd
    # irreversible
    Vmax = parameters["Vmax"]
    S0 = conc[0]
    v = Vmax * S0
    return v
def func_ordered_uni_bi(conc, parameters):
    """
    for alsolase
    """
    Kma = parameters["Kma"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Kiq = parameters["Kiq"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    a = conc[0]
    p = conc[1]
    q = conc[2]
    v = Vmax * (a/Kma)*(1-(p*q/a)/Keq) / (1 + a/Kma + p/Kmp + q/Kmq + a*q/(Kma*Kiq) + p*q/(Kmp*Kmq))
    return v
def func_PYKori(conc, parameters):
    """
    for pyruvate kinase
    """
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]
    pp1 = parameters["pp1"]
    pp2 = parameters["pp2"]
    pp3 = parameters["pp3"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    i1 = conc[4]
    i2 = conc[5]
    i3 = conc[6]
    v = (pp1 + (1 - pp1)/(1 + i1)) * (pp2 + (1 - pp2)/(1 + i2)) * (pp3 + (1 - pp3)/(1 + i3)) * Vmax * (a*b/(Kma*Kmb))/(1.0 + a/Kma+ b/Kmb + a/Kma * b/Kmb)
    return v

def func_PYK(conc, parameters):
    """
    for pyruvate kinase
    """
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]
    Ki1 = parameters["Ki1"]
    Ki2 = parameters["Ki2"]
    Ki3 = parameters["Ki3"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    i1 = conc[4]
    i2 = conc[5]
    i3 = conc[6]
    v = (Ki1/(Ki1+i1))* (Ki2/(Ki2+i2))* (Ki3/(Ki3+i3))* Vmax * (a*b/(Kma*Kmb))/ (1.0 + a/Kma+ b/Kmb + a/Kma * b/Kmb)
    return v

def func_PRK(conc, parameters):
    # Substrate 2
    # Product 2
    # reversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]

    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (Ki/(Ki+ef))* Vmax * (a*b/(Kma*Kmb))/ (1.0 + a/Kma+ b/Kmb + a/Kma * b/Kmb)
    return v


def func_reversible_Michaelis_Menten(conc, parameters):
    # Substrate 1
    # Product 1
    # reversible
    Kmsub = parameters["Kma"]
    Kmpro = parameters["Kmp"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    Sub = conc[0]
    Pro = conc[1]
    v = Vmax * (Sub/Kmsub * (1-(Pro/Sub)/Keq))/(1 + Sub/Kmsub + Pro/Kmpro)
    return v
def func_reversible_Michaelis_Menten_bi_bi(conc, parameters):
    # Substrate 2
    # Product 2
    # reversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    v = Vmax * (a*b/(Kma*Kmb))*(1.0-(p*q/(a*b))/Keq) / ((1.0 + a/Kma + p/Kmp)*(1+ b/Kmb + q/Kmq))
    return v

def func_reversible_Michaelis_Menten_bi_bi_twoisozyme(conc, parameters):
    # Substrate 2
    # Product 2
    # reversible
    Vmax = parameters["Vmax"]
    Kma1 = parameters["Kma1"]
    Kmb1 = parameters["Kmb1"]
    Kmp1 = parameters["Kmp1"]
    Kmq1 = parameters["Kmq1"]
    Vmax1 = parameters["Vmax1"]
    Kma2 = parameters["Kma2"]
    Kmb2 = parameters["Kmb2"]
    Kmp2 = parameters["Kmp2"]
    Kmq2 = parameters["Kmq2"]
    Vmax2 = parameters["Vmax2"]
    Keq = parameters["Keq"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    v = Vmax *(Vmax1 * (a*b/(Kma1*Kmb1))*(1.0-(p*q/(a*b))/Keq) / ((1.0 + a/Kma1 + p/Kmp1)*(1+ b/Kmb1 + q/Kmq1))+\
    Vmax2 * (a*b/(Kma2*Kmb2))*(1.0-(p*q/(a*b))/Keq) / ((1.0 + a/Kma2 + p/Kmp2)*(1+ b/Kmb2 + q/Kmq2)))
    return v

def func_reversible_Michaelis_Menten_bi_uni(conc, parameters):
    # Substrate 2
    # Product 1
    # reversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    v = Vmax * (a*b/(Kma*Kmb))*(1.0-(p/(a*b))/Keq) / ((1.0 + a/Kma + p/Kmp)*(1+ b/Kmb))
    return v

def func_ordered_uni_bi_i(conc, parameters):
    """
    for alsolase
    """
    Kma = parameters["Kma"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Kiq = parameters["Kiq"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    a = conc[0]
    p = conc[1]
    q = conc[2]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))* Vmax * (a/Kma)*(1-(p*q/a)/Keq) / (1 + a/Kma + p/Kmp + q/Kmq + a*q/(Kma*Kiq) + p*q/(Kmp*Kmq))
    return v
def func_PYK_i(conc, parameters):
    """
    for pyruvate kinase
    """
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]
    pp1 = parameters["pp1"]
    pp2 = parameters["pp2"]
    pp3 = parameters["pp3"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    i1 = conc[4]
    i2 = conc[5]
    i3 = conc[6]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))* (pp1 + (1 - pp1)/(1 + i1)) * (pp2 + (1 - pp2)/(1 + i2)) * (pp3 + (1 - pp3)/(1 + i3)) * Vmax * (a*b/(Kma*Kmb))/(1.0 + a/Kma+ b/Kmb + a/Kma * b/Kmb)
    return v

def func_reversible_Michaelis_Menten_i(conc, parameters):
    # Substrate 1
    # Product 1
    # reversible
    Kmsub = parameters["Kma"]
    Kmpro = parameters["Kmp"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    Sub = conc[0]
    Pro = conc[1]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))*  Vmax * (Sub/Kmsub * (1-(Pro/Sub)/Keq))/(1 + Sub/Kmsub + Pro/Kmpro)
    return v
def func_reversible_Michaelis_Menten_bi_bi_i(conc, parameters):
    # Substrate 2
    # Product 2
    # reversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Kmq = parameters["Kmq"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))*  Vmax * (a*b/(Kma*Kmb))*(1.0-(p*q/(a*b))/Keq) / ((1.0 + a/Kma + p/Kmp)*(1+ b/Kmb + q/Kmq))
    return v

def func_reversible_Michaelis_Menten_bi_bi_twoisozyme_i(conc, parameters):
    # Substrate 2
    # Product 2
    # reversible
    Vmax = parameters["Vmax"]
    Kma1 = parameters["Kma1"]
    Kmb1 = parameters["Kmb1"]
    Kmp1 = parameters["Kmp1"]
    Kmq1 = parameters["Kmq1"]
    Vmax1 = parameters["Vmax1"]
    Kma2 = parameters["Kma2"]
    Kmb2 = parameters["Kmb2"]
    Kmp2 = parameters["Kmp2"]
    Kmq2 = parameters["Kmq2"]
    Vmax2 = parameters["Vmax2"]
    Keq = parameters["Keq"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    q = conc[3]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))*  Vmax *(Vmax1 * (a*b/(Kma1*Kmb1))*(1.0-(p*q/(a*b))/Keq) / ((1.0 + a/Kma1 + p/Kmp1)*(1+ b/Kmb1 + q/Kmq1))+\
    Vmax2 * (a*b/(Kma2*Kmb2))*(1.0-(p*q/(a*b))/Keq) / ((1.0 + a/Kma2 + p/Kmp2)*(1+ b/Kmb2 + q/Kmq2)))
    return v

def func_reversible_Michaelis_Menten_bi_uni_i(conc, parameters):
    # Substrate 2
    # Product 1
    # reversible
    Kma = parameters["Kma"]
    Kmb = parameters["Kmb"]
    Kmp = parameters["Kmp"]
    Vmax = parameters["Vmax"]
    Keq = parameters["Keq"]
    a = conc[0]
    b = conc[1]
    p = conc[2]
    Ki = parameters["Ki"]
    ef = conc[-1]
    ppi = parameters["ppi"]
    v = (ppi +(1-ppi)*Ki/(Ki+ef))*  Vmax * (a*b/(Kma*Kmb))*(1.0-(p/(a*b))/Keq) / ((1.0 + a/Kma + p/Kmp)*(1+ b/Kmb))
    return v


#
# Calc functions outside of EnsembleModel class. These functions were required for parallel python
#
def calc_external(conc_list, parameters, reactions, metabolites,reaction_list, metabolite_list, matrix):
    # calc delta of flux levels
    flux = []
    for name in reaction_list:
        conc = [conc_list[metabolites[metname]['id']] for metname in reactions[name]['metabolites']]
        effector = [conc_list[metabolites[metname]['id']]  for metname in reactions[name]['effector']]
        parameters = reactions[name]['parameters']
        func = reactions[name]['func']
        flux.append(func(conc+effector, parameters))
    delta_conc = numpy.dot(matrix, flux)
    #
    for j in range(len(metabolite_list)):
        metabolite = metabolite_list[j]
        #
        # Metabolite with constant concentration
        #
        if metabolites[metabolite]['type'] == 'initial':
            delta_conc[j] = 0.0

    return(delta_conc)

def odeint_external(reactions, metabolites, reaction_list, metabolite_list, matrix, time = 10 , step = 0.1):
    #
    # Solve model outside of the class
    #
    #import mkl
    #mkl.set_num_threads(1)
    initval =  [metabolites[metabolite]['initialconc'] for metabolite in metabolite_list]
    t = numpy.arange(0, time, step)
    result, output =scipy.integrate.odeint(calc_external,initval,t, (reactions, metabolites,reaction_list, metabolite_list, matrix), full_output=1)
    return(result[-10:], output['message'])

def load_metabolic_model_reactions(filename, format = "text"):

    #
    # Distionary for reaction data
    #
    reactions = {}
    #
    # Conter to keep reaction order
    #
    counter = 0
    #
    # Ititialize mode
    #
    mode = "start"

    with open(filename, 'r') as f:
        if format == "text":
            reader = csv.reader(f, delimiter='\t')
        elif format == "csv":
            reader = csv.reader(f, dialect='excel')
        else:
            print("Unknown format!")
            f.close()
            return False

        for i, row in enumerate(reader):
            if "#" in row[0]:
                continue
            row = [item for item in row if "#" not in item]
            if len(row) < 1:
                continue

            if "//" in row[0]:
                if row[0] == "//Reactions":
                    mode = "Reactions"
                    continue
                if row[0] == "//Endreactions":
                    break
                mode = "other"

            if mode == "Reactions":
                if len(row) < 5:
                    continue
                #print row
                rid = row[0].replace(" ", "")
                reaction = row[1].replace(" ", "")
                typer = row[2].replace(" ", "")
                effector = row[3].replace("nothing", "")
                value = row[4]
                reactions[rid] = {
                'reaction' :reaction,
                'type':typer,
                'effector':effector,
                'value':float(value),
                'order':counter
                }
                counter = counter + 1

    return(reactions)

def load_metabolic_model_metabolites(filename, format = "text"):

    metabolites = {}
    counter = 0

    mode = "start"

    with open(filename, 'r') as f:
        if format == "text":
            reader = csv.reader(f, delimiter='\t')
        elif format == "csv":
            reader = csv.reader(f, dialect='excel')
        else:
            print("Unknown format!")
            f.close()
            return False

        for i, row in enumerate(reader):
            #print row
            if "#" in row[0]:
                continue
            row = [item for item in row if "#" not in item]

            if len(row) < 1:
                continue

            if "//" in row[0]:
                if row[0] == "//Metabolites":
                    mode = "Metabolites"
                    continue
                if row[0] == "//Endmetabolites":
                    break
                mode = "other"
            if len(row) < 2:
                continue
            if mode == "Metabolites":
                name = row[0].replace(" ", "")
                initialconc = row[1]
                typem = row[2].replace(" ", "")

                metabolites[name] = {
                'initialconc' :float(initialconc),
                'type':typem,
                'order':counter
                }

                counter = counter + 1
    return(metabolites)

def load_metabolic_model_parameters(filename, format = "text"):

    dic = {}
    counter = 0

    mode = "start"

    with open(filename, 'r') as f:
        if format == "text":
            reader = csv.reader(f, delimiter='\t')
        elif format == "csv":
            reader = csv.reader(f, dialect='excel')
        else:
            print("Unknown format!")
            f.close()
            return False

        for i, row in enumerate(reader):
            #print row

            if len(row) < 1:
                continue

            if "#" in row[0]:
                continue
            row = [item for item in row if "#" not in item]


            if "//" in row[0]:
                if row[0] == "//Parameters":
                    mode = "Parameters"
                    continue
                if row[0] == "//Endparameters":
                    break
                mode = "other"

            if len(row) < 4:
                continue

            if mode == "Parameters":
                #print row
                name = row[0].replace(" ", "")
                parameter = row[1].replace(" ", "")
                lb = row[2]
                value = row[3]
                ub = row[4]
                type = row[5]
                if not name in dic:
                    dic[name] = {}
                dic[name][parameter] = {
                    'lb':float(lb),
                    'ub':float(ub),
                    'value':float(value),
                    'type':type
                }

    return(dic)

def load_metabolic_model(filename, format = "text"):
    """
    Load metabolic model
    """

    parameters = load_metabolic_model_parameters(filename, format)
    reactions = load_metabolic_model_reactions(filename, format)
    metabolites = load_metabolic_model_metabolites(filename, format)

    return(metabolites, reactions, parameters)


class MetabolicModel:
    """
    Class for single metabolic model
    """
    def __init__(self, metabolites, reactions, parameters):
        self.numberofreactions = 0
        self.numberofmetabolites = 0
        self.reactions = {}
        self.metabolites = {}
        self.conc = {}
        self.flux = {}
        self.time = []
        self.metabolite_list = [""] * len(metabolites)
        self.reaction_list = [""] * len(reactions)
        self.parameters = copy.deepcopy(parameters)
        self.initialconc={}
        self.initialflux={}
        #
        # Add metabolites
        #
        #  Set metabolites
        #
        for name in metabolites:
            initialconc = metabolites[name]['initialconc']
            mtype = metabolites[name]['type']
            order = metabolites[name]['order']
            self.add_metabolite(name, initialconc, mtype, order)
            self.metabolite_list[order] = name        # Add reactions
            self.initialconc[name] = metabolites[name]['initialconc']*1.0
        #
        # Set reactions
        #
        for name in reactions:
            value= reactions[name]['value']
            reaction = reactions[name]['reaction']
            rtype = reactions[name]['type']
            order = reactions[name]['order']
            effector = reactions[name]['effector']
            self.add_reaction(name, reaction, rtype , value, effector, order)
            self.initialflux[name] = reactions[name]['value']*1.0
        #
        # Set parameters
        #
        for name in parameters:
            for parameter in parameters[name]:
                self.reactions[name]['parameters'][parameter] = float(parameters[name][parameter]['value'])

        #
        # Construct stoichiometry matrix
        #
        self.construct_stoichiometry_matrix()

    def generate_function(seif, type):
        #
        # Generate functions for rate equations
        #。
        func_tuple = ()
        if type == "Constant":
            reaction_func = func_constant
            paremeters = {'Vmax':1.0}
        func_tuple = func_tuple + (func_constant,)

        if type == "irreversible_Michaelis_Menten":
            reaction_func = func_irreversible_Michaelis_Menten
            paremeters = {'Km':1.0, 'Vmax':1.0}
        func_tuple = func_tuple + (func_irreversible_Michaelis_Menten,)
        if type == "irreversible_Michaelis_Menten_bi_bi":
            reaction_func = func_irreversible_Michaelis_Menten_bi_bi
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kmb':1.0,  'Vmax':1.0}
        func_tuple = func_tuple + (func_irreversible_Michaelis_Menten_bi_bi,)
        if type == "irreversible_Michaelis_Menten_bi_uni":
            reaction_func = func_irreversible_Michaelis_Menten_bi_uni
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmb':1.0,  'Vmax':1.0}
        func_tuple = func_tuple + (func_irreversible_Michaelis_Menten_bi_uni,)


        if type == "irreversible_Michaelis_Menten_i":
            reaction_func = func_irreversible_Michaelis_Menten_i
            paremeters = {'Km':1.0, 'Vmax':1.0, 'Ki':1.0 }
        func_tuple = func_tuple + (func_irreversible_Michaelis_Menten_i,)
        if type == "irreversible_Michaelis_Menten_bi_bi_i":
            reaction_func = func_irreversible_Michaelis_Menten_bi_bi_i
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kmb':1.0,  'Vmax':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_irreversible_Michaelis_Menten_bi_bi_i,)
        if type == "irreversible_Michaelis_Menten_bi_uni_i":
            reaction_func = func_irreversible_Michaelis_Menten_bi_uni_i
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmb':1.0,  'Vmax':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_irreversible_Michaelis_Menten_bi_uni_i,)



        if type == "Mass_action_three_substrate":
            reaction_func = func_mass_action_three_substrate
            paremeters = {'Vmax':1.0}
        func_tuple = func_tuple + (func_mass_action_three_substrate,)
        if type == "Mass_action_two_substrate":
            reaction_func = func_mass_action_two_substrate
            paremeters = {'Vmax':1.0}
        func_tuple = func_tuple + (func_mass_action_two_substrate,)
        if type == "Mass_action_one_substrate":
            reaction_func = func_mass_action_one_substrate
            paremeters = {'Vmax':1.0}
        func_tuple = func_tuple + (func_mass_action_one_substrate,)
        if type == "ordered_uni_bi":
            reaction_func = func_ordered_uni_bi
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kiq':1.0,  'Vmax':1.0, 'Keq':1.0}
        func_tuple = func_tuple + (func_ordered_uni_bi,)

        if type == "PYK":
            reaction_func = func_PYK
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kmb':1.0, 'pp1':1.0, 'pp2':1.0, 'pp3':1.0, 'pp4':1.0, 'qq1':1.0, 'Vmax':1.0}
        func_tuple = func_tuple + (func_PYK,)

        if type == "PRK":
            reaction_func = func_PRK
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kmb':1.0, 'ki1':1.0, 'ki2':1.0, 'ki3':1.0, 'pp4':1.0, 'qq1':1.0, 'Vmax':1.0}
        func_tuple = func_tuple + (func_PRK,)


        if type == "reversible_Michaelis_Menten":
            reaction_func = func_reversible_Michaelis_Menten
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Vmax':1.0, 'Keq':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten,)

        if type == "reversible_Michaelis_Menten_bi_bi_twoisozyme":
            reaction_func = func_reversible_Michaelis_Menten_bi_bi_twoisozyme
            paremeters = {'Kma1':1.0,'Kmp1':1.0, 'Kmq1':1.0, 'Kmb1':1.0,  'Vmax1':1.0, 'Kma2':1.0,'Kmp2':1.0, 'Kmq2':1.0, 'Kmb2':1.0,  'Vmax2':1.0,'Vmax':1.0,  'Keq':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten_bi_bi_twoisozyme,)

        if type == "reversible_Michaelis_Menten_bi_bi":
            reaction_func = func_reversible_Michaelis_Menten_bi_bi
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kmb':1.0,  'Vmax':1.0, 'Keq':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten_bi_bi,)

        if type == "reversible_Michaelis_Menten_bi_uni":
            reaction_func = func_reversible_Michaelis_Menten_bi_uni
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmb':1.0,  'Vmax':1.0, 'Keq':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten_bi_uni,)



        if type == "ordered_uni_bi_i":
            reaction_func = func_ordered_uni_bi_i
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kiq':1.0,  'Vmax':1.0, 'Keq':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_ordered_uni_bi_i,)

        if type == "PYK_i":
            reaction_func = func_PYK_i
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kmb':1.0, 'pp1':1.0, 'pp2':1.0, 'pp3':1.0, 'pp4':1.0, 'qq1':1.0, 'Vmax':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_PYK_i,)

        if type == "reversible_Michaelis_Menten_i":
            reaction_func = func_reversible_Michaelis_Menten_i
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Vmax':1.0, 'Keq':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten_i,)

        if type == "reversible_Michaelis_Menten_bi_bi_twoisozyme_i":
            reaction_func = func_reversible_Michaelis_Menten_bi_bi_twoisozyme_i
            paremeters = {'Kma1':1.0,'Kmp1':1.0, 'Kmq1':1.0, 'Kmb1':1.0,  'Vmax1':1.0, 'Kma2':1.0,'Kmp2':1.0, 'Kmq2':1.0, 'Kmb2':1.0,  'Vmax2':1.0,'Vmax':1.0,  'Keq':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten_bi_bi_twoisozyme_i,)

        if type == "reversible_Michaelis_Menten_bi_bi_i":
            reaction_func = func_reversible_Michaelis_Menten_bi_bi_i
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmq':1.0, 'Kmb':1.0,  'Vmax':1.0, 'Keq':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten_bi_bi_i,)

        if type == "reversible_Michaelis_Menten_bi_uni_i":
            reaction_func = func_reversible_Michaelis_Menten_bi_uni_i
            paremeters = {'Kma':1.0,'Kmp':1.0, 'Kmb':1.0,  'Vmax':1.0, 'Keq':1.0, 'Ki':1.0}
        func_tuple = func_tuple + (func_reversible_Michaelis_Menten_bi_uni_i,)




        if type == 'func_tuple':
            return(func_tuple)
        return(reaction_func, paremeters)

    def add_metabolite(self, name, initialconc, metabolitetype, order):
        """
        Add metabolie to model
        """
        self.metabolites[name] = {'id': order, 'initialconc': initialconc, "type":metabolitetype}
        #self.conc[name] = [initialconc]
        self.conc[name] = []
        self.numberofmetabolites = self.numberofmetabolites + 1
        self.metabolite_list[order] = name

    def add_reaction(self, name, reaction, reactiontype, initialflux, effector, order):
        """
        Add reaction to model。
        """
        self.flux[name] = []
        self.reaction_list[order] = name
        reaction_func, paremeters = self.generate_function(reactiontype)
        self.reactions[name] = {'id': order, 'func': reaction_func,'parameters': paremeters, "type":reactiontype}
        self.reactions[name]['formula'] = reaction
        self.reactions[name]['stoichiometry'] = {}
        self.reactions[name]['initialflux'] = initialflux
        reaction = reaction.replace(" ","")
        reaction = reaction.replace("\t","")
        reaction_temp = re.sub(r'{[0-9\.]+}','', reaction)
        substrate, product = reaction_temp.split("=")
        substrates = substrate.split("+")
        products = product.split("+")
        self.reactions[name]['metabolites'] = substrates + products
        # for stoichiometry matrix
        substrate, product = reaction.split("=")
        substrates = substrate.split("+")
        products = product.split("+")
        for substrate in substrates:
            substrate = substrate.replace("{","")
            substrate_temp = substrate.split("}")
            #print substrate, substrate_temp
            if len(substrate_temp) > 1:
                metabolite = substrate_temp[1]
                self.reactions[name]['stoichiometry'][metabolite] = float(substrate_temp[0]) * -1.0
            else:
                self.reactions[name]['stoichiometry'][substrate] = -1.0
        for product in products:
            product = product.replace("{","")
            product_temp = product.split("}")
            if len(product_temp)  > 1:
                metabolite = product_temp[1]
                self.reactions[name]['stoichiometry'][metabolite] = product_temp[0]
            else:
                self.reactions[name]['stoichiometry'][product] = 1.0
        #
        # Effectors
        #
        effectors = effector.split("+")
        if not effectors[0] == "":
            self.reactions[name]['effector'] = effectors
        else:
            self.reactions[name]['effector'] = []

        self.numberofreactions = self.numberofreactions + 1

    def calc_flux(self, name):
        #
        # Calc flux level for reaction "name" using cerrent metabolic data
        #
        conc = [self.conc[metname][-1] for metname in self.reactions[name]['metabolites']]
        effector = [self.conc[metname][-1] for metname in self.reactions[name]['effector']]
        parameters = self.reactions[name]['parameters']
        func = self.reactions[name]['func']
        return(func(conc+effector, parameters))

    def construct_stoichiometry_matrix(self):
        #
        # construct a stoichiometry matrix
        #
        self.matrix = numpy.zeros((self.numberofmetabolites ,self.numberofreactions))
        for reaction in self.reactions:
            for metabolite in self.reactions[reaction]['stoichiometry']:
                x = self.metabolites[metabolite]['id']
                y = self.reactions[reaction]['id']
                self.matrix[x,y] = self.reactions[reaction]['stoichiometry'][metabolite]
        return(self.matrix)

    def calc_flux_from_conc(self, name, conc_list):
        #
        #Calc flux level for reaction "name" using given metabolic data
        #
        conc = [conc_list[self.metabolites[metname]['id']] for metname in self.reactions[name]['metabolites']]
        effector = [conc_list[self.metabolites[metname]['id']] for metname in self.reactions[name]['effector']]
        parameters = self.reactions[name]['parameters']
        func = self.reactions[name]['func']
        return(func(conc+effector, parameters))

    def calc(self, conc, parameters):
        #
        # calc delta flux from given conc and parameters
        #
        flux = [self.calc_flux_from_conc(reaction, conc) for reaction in self.reaction_list]
        delta_conc = numpy.dot(self.matrix, flux)
        #
        for j in range(self.numberofmetabolites):
            metabolite = self.metabolite_list[j]
            #
            # Metabolite with constance conc
            #
            if self.metabolites[metabolite]['type'] == 'initial':
                delta_conc[j] = 0.0

        return(delta_conc)

    def odeint(self, time = 10, step = 0.1):
        #
        # Solve
        #
        initval =  [self.metabolites[metabolite]['initialconc'] for metabolite in self.metabolite_list]
        t = numpy.arange(0, time, step)
        result, output = scipy.integrate.odeint(self.calc,initval,t, full_output=1)
        length = len(result)

        for i in range(length):
            for j in range(self.numberofmetabolites):
                metabolite = self.metabolite_list[j]
                self.conc[metabolite].append(result[i,j])

            flux = [self.calc_flux(reaction) for reaction in self.reaction_list]
            for j in range(self.numberofreactions):
                reaction = self.reaction_list[j]
                self.flux[reaction].append(flux[j])
        self.time = self.time + list(t)
        return(result, output['message'])

    def odeint_external(self, time = 10, step = 0.1):
        #
        # Solve outside of the model
        #
        result, output = odeint_external(self.reactions, self.metabolites,\
                self.reaction_list, self.metabolite_list, self.matrix,\
                time, step)
        length = len(result)
        for i in range(length):
            for j in range(self.numberofmetabolites):
                metabolite = self.metabolite_list[j]
                self.conc[metabolite].append(result[i,j])

            flux = [self.calc_flux(reaction) for reaction in self.reaction_list]
            for j in range(self.numberofreactions):
                reaction = self.reaction_list[j]
                self.flux[reaction].append(flux[j])
        self.time = self.time  + list(numpy.arange(0.0, time, step))
        return(result, output['message'])




    def set_Vmax_to_given_flux(self, target = []):
        #
        # Set Vmax level for given [S], v and parameteres
        #
        target = set(target)

        for name in self.reactions:
            if len(target) > 0:
                if not name in target:
                    continue
            conc = [self.metabolites[metname]['initialconc'] for metname in self.reactions[name]['metabolites']]
            effector = [self.metabolites[metname]['initialconc'] for metname in self.reactions[name]['effector']]
            parameters = copy.copy(self.reactions[name]['parameters'])
            func = self.reactions[name]['func']
            #
            # Set new Vmax
            #
            if func(conc + effector, parameters) == 0:
                #
                # if reaction rate is zero
                #
                new_Vmax = 0
            else:
                #
                # Calc new Vmax
                #
                new_Vmax = self.reactions[name]['initialflux']/func(conc + effector, parameters) * self.reactions[name]['parameters']['Vmax']
                if new_Vmax >= 0:
                    #
                    # When new_Vmax is OK
                    #
                    pass
                else:
                    #
                    # When new_Vmax is not OK, range of Keq was expanded
                    #
                    multiply = 3.0
                    #print("Range of Keq is expanded in ", name)
                    for i in range(10):
                        parameters['Keq'] = self.reactions[name]['parameters']['Keq'] * multiply
                        new_Vmax = self.reactions[name]['initialflux']/func(conc + effector, parameters) * self.reactions[name]['parameters']['Vmax']
                        #print  name,multiply,new_Vmax,parameters['Keq']
                        if new_Vmax > 0:
                            self.reactions[name]['parameters']['Keq'] = parameters['Keq'] * 1.0
                            self.parameters[name]['Keq']['value'] = parameters['Keq'] * 1.0
                            break
                        else:
                            parameters['Keq'] = self.reactions[name]['parameters']['Keq'] / multiply
                            new_Vmax = self.reactions[name]['initialflux']/func(conc + effector, parameters) * self.reactions[name]['parameters']['Vmax']
                            #print  name,multiply,new_Vmax,parameters['Keq'],"/"
                            if new_Vmax > 0:
                                self.reactions[name]['parameters']['Keq'] = parameters['Keq'] * 1.0
                                self.parameters[name]['Keq']['value'] = parameters['Keq'] * 1.0
                                break
                        multiply = multiply * 3.0
                        #print  name,multiply,new_Vmax
            self.reactions[name]['parameters']['Vmax'] = new_Vmax
            self.parameters[name]['Vmax']['value'] = new_Vmax
            #print "Vmax of ", name, " was ", new_Vmax
            if new_Vmax < 0:
                    print("Vmax of ", name, " was negative value!!!", new_Vmax)


    def set_random_parameter(self, method = "normal"):
        #
        # Set parameters with random value
        #
        for name in self.parameters:
            for parameter in self.parameters[name]:
                lb = self.parameters[name][parameter]['lb']
                ub = self.parameters[name][parameter]['ub']
                type = self.parameters[name][parameter]['type']
                if type == "log":
                    new_value = numpy.random.rand() * (log10(ub)-log10(lb)) + log10(lb)
                    self.reactions[name]['parameters'][parameter] = pow(10,new_value)
                    self.parameters[name][parameter]['value'] = pow(10,new_value)
                elif type == "digital":
                    new_value = int(numpy.random.rand() + 0.5)
                    self.reactions[name]['parameters'][parameter] = new_value
                    self.parameters[name][parameter]['value'] = new_value
                else:
                    new_value = numpy.random.rand() * (ub-lb) + lb
                    self.reactions[name]['parameters'][parameter] = new_value
                    self.parameters[name][parameter]['value'] = new_value

        #
        # Set feasible Keq
        #
        for name in self.reactions:
            conc = [self.metabolites[metname]['initialconc'] for metname in self.reactions[name]['metabolites']]
            effector = [self.metabolites[metname]['initialconc'] for metname in self.reactions[name]['effector']]
            parameters = copy.copy(self.reactions[name]['parameters'])
            func = self.reactions[name]['func']
            newVmax_center = 0
            if func(conc + effector, parameters) != 0:
                newVmax_center = self.reactions[name]['initialflux']/func(conc + effector, parameters)

            if newVmax_center >= 0:
                pass
            else:
                #print "Infeasible Keq in", name
                parameters['Keq'] = self.parameters[name]['Keq']['lb'] * 1.0
                newVmax_lb = 0
                if func(conc + effector, parameters) != 0:
                    newVmax_lb = self.reactions[name]['initialflux']/func(conc + effector, parameters)
                parameters['Keq'] = self.parameters[name]['Keq']['ub'] * 1.0
                newVmax_ub = 0
                if func(conc + effector, parameters) != 0:
                    newVmax_ub = self.reactions[name]['initialflux']/func(conc + effector, parameters)

                lb = self.parameters[name]['Keq']['lb']
                ub = self.parameters[name]['Keq']['ub']
                type = self.parameters[name]['Keq']['type']

                if newVmax_lb > 0:
                    ub = self.parameters[name]['Keq']['value']
                    self.reactions[name]['parameters']['Keq'] = self.parameters[name]['Keq']['lb'] * 1.0
                    self.parameters[name]['Keq']['value'] = self.parameters[name]['Keq']['lb'] * 1.0

                elif newVmax_ub > 0:
                    lb = self.parameters[name]['Keq']['value']
                    self.reactions[name]['parameters']['Keq'] = self.parameters[name]['Keq']['ub'] * 1.0
                    self.parameters[name]['Keq']['value'] = self.parameters[name]['Keq']['ub'] * 1.0

                else:
                    #
                    #　Give up!
                    #
                    #print("Keq value range is infeasible for ", name)
                    continue
                #
                # 1000 times challenge by random value
                #

                for i in range(1000):
                    if type == "log":
                        new_value = numpy.random.rand() * (log10(ub)-log10(lb)) + log10(lb)
                        parameters['Keq']  = pow(10,new_value)
                        newVmax = self.reactions[name]['initialflux']/func(conc + effector, parameters)
                    else:
                        new_value = numpy.random.rand() * (ub-lb) + lb
                        parameters['Keq']  = new_value
                        newVmax = self.reactions[name]['initialflux']/func(conc + effector, parameters)
                    if newVmax > 0:
                        #
                        # if success
                        #
                        self.reactions[name]['parameters']['Keq'] = parameters['Keq'] * 1.0
                        self.parameters[name]['Keq']['value'] = parameters['Keq'] * 1.0



    def activate_reaction(self, reaction, value, parameter= "Vmax", type = "multiply"):
        if type == 'multiply':
            self.reactions[reaction]['parameters'][parameter] = self.reactions[reaction]['parameters'][parameter] * value
        elif type == 'addition':
            self.reactions[reaction]['parameters'][parameter] = self.reactions[reaction]['parameters'][parameter] + value
        elif type == 'set':
            self.reactions[reaction]['parameters'][parameter] = value



    def set_initial_conc(self, metabolite, value, type = "set"):
        if type == 'multiply':
            self.metabolites[metabolite]['initialconc'] = self.metabolites[metabolite]['initialconc'] * value
        elif type == 'addition':
            self.metabolites[metabolite]['initialconc'] = self.metabolites[metabolite]['initialconc'] + value
        elif type == 'set':
            self.metabolites[metabolite]['initialconc'] = value

    def set_initial_flux(self, reaction, value, type = "set"):
        if type == 'multiply':
            self.reactions[reaction]['initialflux'] = self.reactions[reaction]['initialflux'] * value
        elif type == 'addition':
            self.reactions[reaction]['initialflux'] = self.reactions[reaction]['initialflux'] + value
        elif type == 'set':
            self.reactions[reaction]['initialflux'] = value

    def reset_model(self):
        self.time = []
        #
        # Reset parameters
        #
        for name in parameters:
            for parameter in parameters[name]:
                self.reactions[name]['parameters'][parameter] = float(self.parameters[name][parameter]['value'])
        for name in self.reactions:
            #self.flux[name] = [self.reactions[name]['initialflux']]
            self.flux[name] = []
            self.reactions[name]['initialflux'] = self.initialflux[name]
        for name in self.metabolites:
            #self.conc[name] = [self.metabolites[name]['initialconc']]
            self.metabolites[name]['initialconc'] = self.initialconc[name]
            self.conc[name] = []

    def update_model(self):
        #
        # Update model by the current state
        #
        for name in parameters:
            for parameter in parameters[name]:
                self.parameters[name][parameter]['value'] = float(self.reactions[name]['parameters'][parameter])
        for name in self.reactions:
            self.initialflux[name] = self.reactions[name]['initialflux']
        for name in self.metabolites:
            self.initialconc[name] = self.metabolites[name]['initialconc']

    def update_initials(self):
        #
        # Update initial values of metabolic model by the current state
        #
        for name in self.reactions:
            self.reactions[name]['initialflux'] = self.flux[name][-1]
        for name in self.metabolites:
            self.metabolites[name]['initialconc'] = self.conc[name][-1]


    def clear_result(self):
        #
        self.time = []
        #
        # Clear dvelopment history
        #
        for name in self.reactions:
            self.flux[name] = []
        for name in self.metabolites:
            self.conc[name] = []

    def check_convergence(self, thres = 0.01):
        #
        #
        #
        non_passed = {}
        for name in self.metabolites:
            if self.conc[name][-1] > 0.1:
                if self.conc[name][-2]/self.conc[name][-1] < (1 - thres):
                    non_passed[name] = 1;
                if self.conc[name][-2]/self.conc[name][-1] > (1 + thres):
                    non_passed[name] = 1;
            else:
                if abs(self.conc[name][-2] - self.conc[name][-1]) > (0.001):
                    non_passed[name] = 1;
        return([name for name in non_passed])



    def show_parameter(self):
        #
        for name in sorted(self.reactions):
            for parameter in sorted(self.reactions[name]['parameters']):
                print(name+"\t"+parameter+"\t"+str(self.reactions[name]['parameters'][parameter]))
    def show_flux(self, reaction):
        for i in range(len(self.time)):
            print(str(self.time[i])+"\t"+str(reaction)+"\t"+str(self.flux[reaction][i]))
    def show_fluxes(self):
        for name in self.reaction_list:
            print(str(name)+"\t"+str(self.flux[name][-1]))
    def show_conc(self, metabolite):
        for i in range(len(self.time)):
            print(str(self.time[i])+"\t"+str(metabolite)+"\t"+str(self.conc[metabolite][i]))
    def show_concs(self):
        for name in self.metabolite_list:
            print(str(name)+"\t"+str(self.conc[name][-1]))

    def get_parameters(self, time, step):
        return((self.reactions, self.metabolites, self.reaction_list, self.metabolite_list, self.matrix,\
            time, step))

    def copy_parameter(self, reactionfrom, parameterfrom,reactionto, parameterto):

        self.reactions[reactionto]['parameters'][parameterto] = self.reactions[reactionfrom]['parameters'][parameterfrom]
        self.parameters[reactionto][parameterto]['value'] = self.reactions[reactionfrom]['parameters'][parameterfrom]


    def delta(self, type = "view"):
        conc =  [self.conc[metabolite][-1] for metabolite in self.metabolite_list]
        flux = [self.calc_flux_from_conc(reaction, conc) for reaction in self.reaction_list]
        delta_conc = numpy.dot(self.matrix, flux)
        #
        #
        SS = 0
        for j in range(self.numberofmetabolites):
            metabolite = self.metabolite_list[j]
            if self.metabolites[metabolite]['type'] == 'free':
                SS = SS + delta_conc[j]*delta_conc[j]
                if type == "view":
                    print(metabolite, delta_conc[j])
        return(SS**0.5)

    def lsa(self, type = "calc"):
        """
        Linar stability analysis
        """
        conc =  [self.conc[metabolite][-1] for metabolite in self.metabolite_list]
        flux = [self.calc_flux_from_conc(reaction, conc) for reaction in self.reaction_list]
        delta_conc = numpy.dot(self.matrix, flux)
        #
        #
        delta = 0.0001
        dxdu = numpy.zeros((len(delta_conc), len(delta_conc)))
        for j in range(len(delta_conc)):
            conc_temp = list(conc)
            conc_temp[j] = conc_temp[j] + delta
            flux_temp = [self.calc_flux_from_conc(reaction, conc_temp) for reaction in self.reaction_list]
            delta_conc_temp = numpy.dot(self.matrix, flux_temp)
            dxdu[:,j] = (delta_conc_temp - delta_conc)/delta
        # Calc eigen values
        la, v = numpy.linalg.eig(dxdu)
        if type=="view":
            print(la)
        check = 0
        for i in la:
            if i.real > 0.000001:
                check = 1

        return(check)


    def save_timecourse(self, filename, format = "text"):
        """
        Save MDV data to the text/csv file

        Parameters

        """
        #
        # preparation of data
        #
        Data = []
        string = ["time"]
        for i in range(len(self.metabolites)):
            metabolite = self.metabolite_list[i]
            string.append(metabolite)
        for i in range(len(self.reactions)):
            reaction = self.reaction_list[i]
            string.append(reaction)
        Data.append(string)
        factor = 1
        if len(self.time) > 3000:
            factor = int(len(self.time)/1000)
            print("Factor is", factor)
        for j in range(len(self.time)):
            if j % factor != 0:
                continue
            string = [self.time[j]]
            for i in range(len(self.metabolites)):
                metabolite = self.metabolite_list[i]
                string.append(self.conc[metabolite][j])
            for i in range(len(self.reactions)):
                reaction = self.reaction_list[i]
                string.append(self.flux[reaction][j])
            Data.append(string)

        with open(filename, 'w') as f: #fo rpython 3
            import csv
            if format == "text":
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(Data)
            if format == "csv":
                writer = csv.writer(f, dialect='excel')
                writer.writerows(Data)
        f.close()
        return True
    def add_allosteric_regulation(self, reaction, effector, mode = "i"):
        """
        Save MDV data to the text/csv file

        Parameters

        """
        self.reactions[reaction]["effector"].append(effector)
        reactiontype =  self.reactions[reaction]['type']+"_i"
        reaction_func, paremeters = self.generate_function(reactiontype)
        self.reactions[reaction]['func'] = reaction_func
        #self.reactions[reaction]['parameters'] = paremeters
        self.reactions[reaction]['type'] = reactiontype
        return True

class EnsembleModel:
    def __init__(self, metabolites, reactions, parameters, number, numberofcpus=3):
        self.models = []
        self.model_use = []
        self.model_output = []
        self.numberofcpus = numberofcpus
        for i in range(number):
            self.models.append(MetabolicModel(metabolites, reactions, parameters))
            self.model_use.append('use')
            self.model_output.append('new')

    def model_selection(self, type, parameter, direction, thres):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                if type == 'conc':
                    value = emodel.models[i].conc[parameter][-1]
                elif type == 'flux':
                    value = emodel.models[i].flux[parameter][-1]
                if direction == 'morethan':
                    if value <= thres:
                        self.model_use[i] = 'notuse'
                        print("model ", i, " was removed due to ",parameter, value, "was less than", thres)
                elif direction == 'lessthan':
                    if value >= thres:
                        self.model_use[i] = 'notuse'
                        print("model ", i, " was removed due to ",parameter, value, "was more than", thres)

    def model_selection_by_finish_condition(self):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                if self.model_output[i] != 'Integration successful.':
                    self.model_use[i] = 'notuse'
                    print("model ", i, " was removed due to ", self.model_output[i])
                #
                # Check metabolite conc
                #
                for metabolite in self.models[i].metabolite_list:
                    if self.models[i].conc[metabolite][-1] < 0:
                        self.model_use[i] = 'notuse'
                        print("model ", i, " was removed due to negative ", metabolite, " level: " ,self.models[i].conc[metabolite][-1])
                        break


    def set_Vmax_to_given_flux(self, target = []):
        for i in range(len(self.models)):
            self.models[i].set_Vmax_to_given_flux(target)

    def odeint(self, time = 10, step = 0.01):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                result, output = self.models[i].odeint(time, step)
                self.model_output[i] = output
    def odeint_external(self, time = 10, step = 0.01):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                result, output = self.models[i].odeint_external(time, step)
                self.model_output[i] = output

    def check_convergence(self, thres = 0.01, type = "delta"):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                if type == "delta":
                    result = self.models[i].delta(type= "no")
                    if result > thres:
                        self.model_use[i] = "notuse"
                        print("model ", i, " was removed due to inconvergence of ", result)
                else:
                    result = self.models[i].check_convergence( thres = 0.01 )
                    if len(result) > 0:
                        self.model_use[i] = "notuse"
                        print("model ", i, " was removed due to inconvergence of ", result)

    def check_stability(self, type = "calc"):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                    result = self.models[i].lsa(type)
                    if result == 1:
                        self.model_use[i] = "notuse"
                        print("model ", i, " was removed due to instability")


    def reset_model(self):
        for i in range(len(self.models)):
            self.models[i].reset_model()
            self.model_output[i] = "new"
            self.model_use[i] = "use"


    def clear_result(self):
        for i in range(len(self.models)):
            self.models[i].clear_result()

    def update_initials(self):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                self.models[i].update_initials()

    def update_models(self):
        for i in range(len(self.models)):
            self.models[i].update_model()

    def set_random_parameter(self):
        for i in range(len(self.models)):
            self.models[i].set_random_parameter()
    def activate_reaction(self, reaction, times, parameter = "Vmax", type = "multiply"):
        for i in range(len(self.models)):
            if type == "set":
                self.models[i].activate_reaction(reaction, times, parameter, "set")
            elif type == "addition":
                self.models[i].activate_reaction(reaction, times, parameter, "addition")
            else:
                self.models[i].activate_reaction(reaction, times, parameter, "multiply")

    def set_initial_flux(self, reaction, times, type = "set"):
        for i in range(len(self.models)):
            self.models[i].set_initial_flux(reaction, times, type = "set")

    def set_initial_conc(self, metabolite, times, type = "set"):
        for i in range(len(self.models)):
            self.models[i].set_initial_conc(metabolite, times, type = "set")

    def show_flux(self, reaction):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                print(str(i)+"\t"+str(reaction)+"\t"+str(self.models[i].flux[reaction][-1]))
    def show_conc(self, metabolite):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                print(str(i)+"\t"+str(metabolite)+"\t"+str(self.models[i].conc[metabolite][-1]))
    def show_delta(self):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                print(str(i)+"\tDelta\t"+str(self.models[i].delta(type = "no")))

    def show_parameter(self, reaction, parameter):
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                print(str(i)+"\t"+str(self.models[i].reactions[reaction]['parameters'][parameter]))
    def odeint_parallel(self, time = 10, step = 0.01):
        t = numpy.arange(0, time, step)

        job_server = pp.Server(ncpus = self.numberofcpus, ppservers=())
        print(job_server.get_active_nodes())
        print("Starting pp with", job_server.get_ncpus(), "workers")
        jobs = []
        #
        # tuples of fuctions to send pp
        #
        external_functions = (calc_external,) + self.models[0].generate_function(type = "func_tuple")
        for i in range(len(self.models)):
            if self.model_use[i] == "use":
                parameters = self.models[i].get_parameters(time, step)
                jobs.append([i, job_server.submit(odeint_external, parameters, external_functions,
                ("numpy","scipy.integrate"))])
                #print("append", i)
        for j, job in jobs:
            results = job()
            result = results[0]
            output = results[1]
            self.model_output[j] = output
            #print(j,output)

            length = len(result)
            for i in range(length):
                for k in range(self.models[j].numberofmetabolites):
                    metabolite = self.models[j].metabolite_list[k]
                    self.models[j].conc[metabolite].append(result[i,k])
                flux = [self.models[j].calc_flux(reaction) for reaction in self.models[j].reaction_list]
                for k in range(self.models[j].numberofreactions):
                    #Retrive reaction name
                    reaction = self.models[j].reaction_list[k]
                    #Retrive reaction value
                    self.models[j].flux[reaction].append(flux[k])
            self.models[j].time = self.models[j].time + list(t)
        job_server.print_stats()
        job_server.destroy()
        return()



    def save_timecourse(self, filename, format = "text"):
        """
        Save time course data to the text/csv file

        """
        #
        # preparation of data
        #
        Data = []
        string = ["time"]
        for model_number in range(len(self.models)):
            for i in range(len(self.models[model_number].metabolites)):
                metabolite = self.models[model_number].metabolite_list[i]
                string.append(str(model_number)+"_"+metabolite)
            for i in range(len(self.models[model_number].reactions)):
                reaction = self.models[model_number].reaction_list[i]
                string.append(str(model_number)+"_"+reaction)
        Data.append(string)
        factor = 1
        if len(self.models[0].time) > 3000:
            factor = int(len(self.models[0].time)/1000)
            print("Factor is", factor)
        for j in range(len(self.models[0].time)):
            if j % factor != 0:
                continue
            string = [self.models[0].time[j]]
            for model_number in range(len(self.models)):
                for i in range(len(self.models[model_number].metabolites)):
                    metabolite = self.models[model_number].metabolite_list[i]
                    string.append(self.models[model_number].conc[metabolite][j])
                for i in range(len(self.models[model_number].reactions)):
                    reaction = self.models[model_number].reaction_list[i]
                    string.append(self.models[model_number].flux[reaction][j])
            Data.append(string)

        with open(filename, 'w') as f:
            import csv
            if format == "text":
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(Data)
            if format == "csv":
                writer = csv.writer(f, dialect='excel')
                writer.writerows(Data)
        f.close()
        return True
    def save_results(self, filename, key, format = "text"):
        """
        Save Parameter data to the text/csv file

        """
        #
        # preparation of data
        #
        reactions = [reaction for reaction in self.models[0].reactions]
        metabolites = [metabolite for metabolite in self.models[0].metabolites]
        string = []
        for reaction in sorted(reactions):
            for i in range(len(self.models)):
                if self.model_use[i] == "use":
                    if len(self.models[i].flux[reaction]) > 0:
                        string.append([str(key), "reaction", reaction, "flux", str(i), str(self.models[i].flux[reaction][-1])])
        for metabolite in sorted(metabolites):
            for i in range(len(self.models)):
                if self.model_use[i] == "use":
                    if len(self.models[i].conc[metabolite]) > 0:
                        string.append([str(key), "met", metabolite, "conc", str(i), str(self.models[i].conc[metabolite][-1])])
        for reaction in sorted(reactions):
            for i in range(len(self.models)):
                if self.model_use[i] == "use":
                    for parameter in sorted(self.models[i].reactions[reaction]['parameters']):
                        string.append([str(key), "parameter", reaction, parameter, str(i), str(self.models[i].reactions[reaction]['parameters'][parameter])])
        for reaction in sorted(reactions):
            for i in range(len(self.models)):
                if self.model_use[i] == "use":
                    string.append([str(key), "initialflux", reaction, "flux", str(i), str(self.models[i].reactions[reaction]['initialflux'])])
        for metabolite in sorted(metabolites):
            for i in range(len(self.models)):
                if self.model_use[i] == "use":
                    string.append([str(key), "initialconc", metabolite, "conc", str(i), str(self.models[i].metabolites[metabolite]['initialconc'])])
        for reaction in sorted(reactions):
            for i in range(len(self.models)):
                if self.model_use[i] == "use":
                    for parameter in self.models[i].parameters[reaction]:
                        string.append([str(key), "originalparameter", reaction, parameter, str(i), str(self.models[i].parameters[reaction][parameter]['value'])])
        with open(filename, 'w') as f:#for Python 3
            if format == "text":
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(string)
            if format == "csv":
                writer = csv.writer(f, dialect='excel')
                writer.writerows(string)
        f.close()
        return True
    def load_parameters(self, filename, format = "text"):
        """
        Load Parameter data from the text/csv file
        """
        for i in range(len(self.models)):
            self.model_use[i] = 'notuse'

        with open(filename, 'r') as f:
            if format == "text":
                reader = csv.reader(f, delimiter='\t')
            elif format == "csv":
                reader = csv.reader(f, dialect='excel')
            else:
                print("Unknown format!")
                f.close()
                return False

            for i, row in enumerate(reader):
                row = [item for item in row]

                if len(row) !=6:
                    continue

                if row[1] == "initialconc":
                    metabolite = row[2]
                    model = int(row[4])
                    value = float(row[5])
                if row[1] == "initialflux":
                    reaction = row[2]
                    model = int(row[4])
                    value = float(row[5])
                if row[1] == "originalparameter":
                    reaction = row[2]
                    parameter = row[3]
                    model = int(row[4])
                    value = float(row[5])
                    self.models[model].parameters[reaction][parameter]['value'] = value
                    self.model_use[model] = 'use'
        f.close()
        return True

    def truncate(self):
        #
        # Remove "notuse" model
        #
        usingmodel = [i for i in range(len(self.models)) if self.model_use[i] == "use"]
        self.model_use = ["use" for i in usingmodel]
        self.models = [self.models[i] for i in usingmodel]
    def report_number_of_models(self):
        #
        # Report number of model in the current ensemble
        #
        usingmodel = [i for i in range(len(self.models)) if self.model_use[i] == "use"]
        print("Current ensemble contains", len(usingmodel), " models.")
        return(len(usingmodel))
    def copy_parameter(self, reactionfrom, parameterfrom,reactionto, parameterto):
        #
        # Copy parameter
        #
        for i in range(len(self.models)):
            self.models[i].copy_parameter(reactionfrom, parameterfrom,reactionto, parameterto)
    def add_allosteric_regulation(self, reaction, effector, mode = "i"):
        """
        Save MDV data to the text/csv file

        Parameters

        """
        for i in range(len(self.models)):
            self.models[i].add_allosteric_regulation(reaction, effector, mode = "i")

if __name__ == '__main__':
    # Read metabolic model definition file
    #
    metabolites, reactions, parameters = load_metabolic_model('PCC6803_200116.txt')
    #
    # 酵素活性バージョン
    #
    import copy
    import itertools
    import statistics
    cpus = 23
    parallel = 1000 # Number of metabolic models in one batch

    #
    # [S]Et-g data
    #
    metabolites_dic = {"G6P": 0.484479487,"F6P": 0.279082024,"FBP": 0.41420439,"DHAP": 0.30185355,"3PG": 25.86966222,"2PG": 10.71674593,
    "PEP": 12.78475698,"6PG": 0.047958882,"Ru5P": 0.01797026,"R5P": 0.020772371,"S7P": 0.030301644,"RuBP": 1.295477693,
    "AcCoA": 0.903188205,"citrate": 2.883932814,"isoCitrate": 3.52005633,
    "Succinate": 0.162060329,"Malate": 0.302218554,"NADP": 2.074724629,"NAD": 0.153243178,"ATP": 1.052973628}

    reaction_list=['prk','gpm','pgk','pyk','pdh','zwf','gnd','rpe','rpi','tkl1','tkl2','tal','aldo','glp','rbc','glt','can','icd','sdh','fum','mdh','ppc','mae','gabD','icl','mls','fbp','pgi','pfk','fba','tpi','gap','eno']
    effector_list=['nothing','citrate','isoCitrate','BPG','E4P','X5P','SBP','G6P','F6P','FBP','DHAP','G3P','3PG','2PG','PEP','PYR','6PG','Ru5P','R5P','S7P','RuBP','AcCoA','2KG','Succinate','Fumarate','Malate','OAA','NADP','NADPH','NAD','NADH','ADP','ATP']

    for (reaction, effector) in itertools.product(reaction_list,effector_list):
        for k in range(0,10):
            print(reaction, effector, k)
            filename ="goodmodels"+str(k)+".txt"
            emodel = EnsembleModel(metabolites,reactions, parameters, parallel, numberofcpus = cpus)
            # Load parameters
            emodel.load_parameters(filename)
            # Reset
            emodel.reset_model()
            #
            # Add allosteric regulation
            #
            if not effector == 'nothing':
                emodel.add_allosteric_regulation(reaction, effector, "i")
                #
                # Modify Vmax to keep initial flux level
                #
                emodel.set_Vmax_to_given_flux(reaction)
            else:
                nothinglist = list(range(1000))
            #
            # [E]auto => [E]Et-g
            #
            emodel.activate_reaction('GlycoDegra',0, type="set")
            emodel.activate_reaction('pdc_constant',0.1127349177992005, type="set")
            emodel.activate_reaction('biomass_auto',0.0129627121234772, type="set")
            #[EEt-g]/[Eauto]
            emodel.activate_reaction('PS2',0.317522436674765)
            emodel.activate_reaction('ATPsynthase',0.640095808432722)
            emodel.activate_reaction('PS1',0.479520247860281)
            emodel.activate_reaction('Cyclic1',0.688279597579098)
            emodel.activate_reaction('Cyclic2',0.688279597579098)

            emodel.activate_reaction('pgi',0.6387550187002)
            emodel.activate_reaction('fbp',0.34833832603842)
            emodel.activate_reaction('pfk',0.452058555153743)
            emodel.activate_reaction('fba',1.07752454705741)
            emodel.activate_reaction('tpi',0.716668406900194)
            emodel.activate_reaction('gap',0.595324281205730, parameter = 'Vmax1')
            emodel.activate_reaction('gap',0.606370809964625, parameter = 'Vmax2')
            emodel.activate_reaction('pgk',0.929872087223499)
            emodel.activate_reaction('gpm',0.781539456006131)
            emodel.activate_reaction('eno',0.765863656898541)
            emodel.activate_reaction('pyk',0.769655896989666)
            emodel.activate_reaction('pdh',0.701521558438302)
            emodel.activate_reaction('zwf',0.606895793756199)
            emodel.activate_reaction('gnd',0.443462187810197)
            emodel.activate_reaction('rpe',1.13007656266981)
            emodel.activate_reaction('rpi',0.719638949797989)
            emodel.activate_reaction('tkl1',0.805515127926094)
            emodel.activate_reaction('tkl2',0.805515127926094)
            emodel.activate_reaction('tal',0.684683951133748)
            emodel.activate_reaction('glp',0.717791920316889)
            emodel.activate_reaction('prk',0.68397405017847)
            emodel.activate_reaction('rbc',0.880701769835223)
            emodel.activate_reaction('icd',0.743110948833053)
            emodel.activate_reaction('sdh',0.544201529672118)
            emodel.activate_reaction('fum',0.609019252373934)
            emodel.activate_reaction('ppc',0.48902208509086)
            emodel.activate_reaction('mae',0.480778036654116)
            emodel.activate_reaction('gabD',0.32668353005897)

            #
            # Solve
            #
            emodel.odeint_parallel(time = 1000, step = 10)
            #
            # model selection
            #
            emodel.model_selection_by_finish_condition()
            emodel.check_convergence(thres = 0.01)
            remains = emodel.report_number_of_models()
            if remains == 0:
                continue
            #
            emodel.update_initials()
            emodel.clear_result()
            #
            # Solve again
            #
            emodel.odeint_parallel(time = 100000, step = 10)
            #
            # model selection
            #
            emodel.model_selection_by_finish_condition()
            emodel.check_convergence(thres = 0.00001)
            emodel.check_stability()
            remains = emodel.report_number_of_models()
            #
            # Evaluation by RSS
            #
            rss_list = []
            rss_log = ""
            rss_counter = 0
            for j in range(len(emodel.models)):

                if emodel.model_use[j] == "use":
                    rss = 0
                    for metabolite in metabolites_dic:
                        if len(emodel.models[j].conc[metabolite]) > 0:
                            conc = emodel.models[j].conc[metabolite][-1]
                            rss =  rss + (conc - metabolites_dic[metabolite]) * (conc - metabolites_dic[metabolite])
                    rss_list.append(rss)
                    rss_log = rss_log + reaction+"\t"+effector+"\t"+str(j)+"\t"+str(rss)+"\n"
                    if effector == "nothing":
                        nothinglist[j] = rss
                    if nothinglist[j] > rss:
                        rss_counter = rss_counter + 1
            if effector == "nothing":
                print(nothinglist)
            print(reaction, effector, remains)
            if remains == 0:
                continue
            #
            # preparation of data
            #
            with open("log_ar.txt", "a") as file:
                rssmean = "nd"
                if remains > 0:
                    rssmean = sum(rss_list) / len(rss_list)
                if len(rss_list) == 0:
                    med = 10000
                else:
                    med = statistics.median(rss_list)
                text = str(k)+"\t"+reaction+"\t"+effector+"\t"+str(len(rss_list))+"\t"+str(rssmean)+"\t"+str(rss_counter)+"\t"+str(med)+"\n"
                file.write(text)
            with open("rsslog_ar.txt", "a") as file:
                file.write(rss_log)


