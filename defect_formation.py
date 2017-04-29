
# coding: utf-8

# In[695]:

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib.ticker import MultipleLocator
import pandas as pd
get_ipython().magic('matplotlib inline')


# In[413]:

# unit eV
# computation restuls
# specify the VBM and CBM
# also VBM value in VASP output
#
VBM = 0
CBM = 3.754
absolute_VBM = 3.014
# polaron self-trapping energy
elec_trap = 0.75
hole_trap = 0.86


# In[428]:

#get_defect_info('V_Mn^')


# In[781]:
#
# set up a dictionary to defect_type, total_energy, correction
# 'Ni_Mn^+':[-1156.15,0.15,1]
# total energy: outputs of VASP
# correction term: Freysoldt correction (Phys. Rev. Lett. 102, 016402)
# The way to automaticly obtain it is under development
#
# format: defect_name:[total energy, correction, charge state]
#
data_Energy = {'bulk':[-1161.146,0.,0],'Na_Mn^0':[-1146.095,0.,0],'Mg_Na-V_Na^0':[-1159.163,0.0,0],'V_Na^-':[-1152.645,0.234,-1],'V_Na^0':[-1156.714,0.0,0],'Mn_Na^0':[-1171.743,0.,0],'Mn_Na^+':[-1177.607,0.185,1],'Mn_Na^2+':[-1181.664,0.697,2],'Mn_Na-Na_Mn^0':[-1158.227,0.,0],'Mn_Na-Na_Mn^-':[-1153.676,0.187,-1],'V_Mn^3-':[-1127.904,1.879,-3],'V_Mn^2-':[-1132.523,0.984,-2],'V_Mn^-':[-1137.128,0.355,-1],'V_Mn^0':[-1140.706,0.,0],'H_i^0':[-1163.966,0.,0],'H_i^+':[-1169.478,0.22,1],'V_O^2+':[-1159.304,0.578,2],'V_O^+':[-1154.655,0.076,1],'V_O^0':[-1149.931,0.,0]}


# In[807]:

data_Energy_DF=pd.DataFrame.from_dict(data_Energy,orient='index')
#obt_TolEnergy_correction('Na_Mn^2-',[-1137.182,0.776,-2],data_Energy)
data_Energy_DF.rename(columns={0:'Total energy (eV)',1:'Correction (eV)',2:'Charge state'}, inplace=True)
#data_Energy_DF


# In[808]:

#data_Energy_DF.reset_index()


# In[809]:

#data_Energy


# In[493]:

# Formation enthalpys
Formation_enthalpy = {'NaMnO2':-7.8395,'GeO2':4.988,'H2O':-2.6655,'TiO2':-9.175333333,'MgO':-5.714,'Fe2O3':-8.8815,'NiO':-2.978,'CoO':-3.266,'Al2O3':-16.072}


# In[3]:

# Na rich
ChemPot1 = {'Na':-1.4975,'Mn':-13.1725,'O':-10.8065,'H':-3.9125,'Ge':-5.353,'Ti':-10.997,'Mg':-3.4555,'Fe':-11.5785,'Ni':-6.36375,'Co':-9.0685,'Al':-6.441}
# Na poor
ChemPot2 = {'Na':-4.0695,'Mn':-16.7795,'O':-7.7165,'H':-4.89175,'Ge':-8.927333333,'Ti':-17.177,'Mg':-6.5455,'Fe':-14.95875,'Ni':-8.63475,'Co':-11.6275,'Al':-11.076} 


# In[811]:

comp_formation_energy(ChemPot1,'V_Mn^3-')


# In[505]:

# about functions to do data analysis
#
# add_TolEnergy_correction_state(str defect_name, list [tot energy, correction, charge state])
# return type: None; just add to data_Energy dictionary
#
# get_defect_info(str defect_name)
# return type: dictionary {elem1:num, elem2:num,...}. This will relate to the formation energy calculation.
#
# comp_formation_energy(dict chem_pot, str defect_name)
# return a float number: the formation energy of the defect.
#
# comp_energy_group(dict chem_pot, str element_char)
# return a dictionary {str defect_name: float formation energy}
# to calculate the formation energy of a group of defects containing the same element
#
# comp_energy_all(dict chem_pot)
# return a dictionary {str defect_name: float formation energy}
# for all the defects
#
# ObtChargeState(str defect_name)
# return the charge state of the defect


# In[503]:

# this function can add the newly obtained results to the data_energy dictionary
def add_TolEnergy_correction_state(defect_str,energy_correction_state_list):
    global data_Energy
    if defect_str in data_store_dic:
        data_Energy[defect_str][0]=energy_correction_list[0]
        data_Energy[defect_str][1]=energy_correction_list[1]
        data_Energy[defect_str][2]=energy_correction_list[2]
        print('Information about '+defect_str+' updated!')
    else:
        data_Energy[defect_str]=energy_correction_list
        print(defect_str+' added!')


# In[443]:

def get_defect_info(defect_name):
    # not include charge state information
    return_dic = {}
    # for example return_dic = {'Na':1} for V_Na
    name_pure = defect_name.split('^')[0]
    if '-' in name_pure:
        seperated = name_pure.split('-')
    else:
        seperated = [name_pure]
    for comp_defect in seperated:
        part1 = comp_defect.split('_')[0]
        part2 = comp_defect.split('_')[1]
        if ord(part2[0]) >= ord('0') and ord(part2[0]) <= ord('9'):
            num_part2 = ord(part2[0])-ord('0')
            real_part2 = part2[1:]
        else:
            num_part2 = 1
            real_part2 = part2
        if part1 != 'V':
            if part1 not in return_dic :
                return_dic[part1] = -1
            else:
                return_dic[part1] -= 1
        if real_part2 !='i':
            if real_part2 not in return_dic:
                return_dic[real_part2] = num_part2
            else:
                return_dic[real_part2] += 1 
    return return_dic


# In[504]:

#get_defect_info('V_Na-Mn_Na^0')


# In[489]:
#
# compute the formation energy of a specific defect under certain chem_pot
# input: chemical potential, defect name
# input_type: dictionary, string
# output: formation energy
# output_type: float
#
def comp_formation_energy(chem_pot, defect_name):
    global data_Energy
    if defect_name not in data_Energy:
        print('Total energy of the defect does not exist!')
        return
    defect_info = get_defect_info(defect_name)
    for str_element in defect_info:
        if str_element not in chem_pot:
            print('The chemical potential of '+str_element+' is not available!')
            return
    form_energy = data_Energy[defect_name][0]-data_Energy['bulk'][0]+data_Energy[defect_name][1]+data_Energy[defect_name][2]*absolute_VBM
    for str_element in defect_info:
        form_energy = form_energy + defect_info[str_element]*chem_pot[str_element]
    return form_energy


# In[509]:

# Output format:
# Formation_Mg_1={'Mg_Na^0':1.245,'Mg_Mn^0':1.625}
#
def comp_energy_group(chem_pot, element_char):
    group_form_energy = {}
    for defect_key in data_Energy:
        if element_char in defect_key:
            if defect_key not in group_form_energy:
                group_form_energy[defect_key]=comp_formation_energy(chem_pot,defect_key)
            else:
                group_form_energy[defect_key]=comp_formation_energy(chem_pot,defect_key)
                print('Duplicated defects:', defect_key)
    return group_form_energy


# In[510]:

comp_energy_group(ChemPot1,'Mn')


# In[511]:
####################################################################
# compute formation energies of different types of defects
# input: the choice of chemical potential
# input_type: dictionary
# output: formation energy of point defects
# output_type: dictionary
####################################################################
def comp_energy_all(chem_pot):
    all_form_energy = {}
    global data_Energy
    for defect_key in data_Energy:
        #if defect_key not in all_form_energy:
        if defect_key != 'bulk':
            all_form_energy[defect_key]=comp_formation_energy(chem_pot,defect_key)
        #else:
        #all_form_energy[defect_key]=comp_formation_energy(data_store_dic,chem_pot,defect_key)
        #print('Duplicated defects:', defect_key)
    return all_form_energy


# In[785]:
#
# convert the dictionary of formation energy to a DataFrame
#
form_Energy_DF=pd.DataFrame.from_dict(comp_energy_all(ChemPot1),orient='index')
form_Energy_DF2=pd.DataFrame.from_dict(comp_energy_all(ChemPot2),orient='index')
form_Energy_DF.rename(columns={0:'Formation enegy (eV) (Na-rich)'},inplace=True)
form_Energy_DF2.rename(columns={0:'Formation enegy (eV) (Na-poor)'},inplace=True)
form_Energy_DF_total=pd.concat([form_Energy_DF, form_Energy_DF2], axis=1)
form_Energy_DF_total.reset_index(inplace=True)
form_Energy_DF_total.sort_values(by='index',axis=0,inplace=True)
form_Energy_DF_total['Defect type'] = form_Energy_DF_total['index'].map(lambda x: x.split('^')[0])
form_Energy_DF_total['Charge state'] = form_Energy_DF_total['index'].map(lambda x: x.split('^')[1])
form_Energy_DF_total.set_index(['Defect type','index'],inplace=True)
form_Energy_DF_total.reset_index(inplace=True)
form_Energy_DF_total.rename(columns={'Defect type': 'Defect','index':'Defect type'},inplace=True)
#form_Energy_DF_total.groupby('Defect')
arrays = [[form_Energy_DF_total['Defect'][i] for i in range(31)],[form_Energy_DF_total['Defect type'][i] for i in range(31)]]
index = pd.MultiIndex.from_arrays(arrays, names=['Defect', 'Defect type'])
# ,'Formation enegy (eV) (Na-rich)','Charge state']
#form_Energy_DF_total['Formation enegy (eV) (Na-rich)']
array1 = [form_Energy_DF_total['Formation enegy (eV) (Na-poor)'][i] for i in range(31)]
array2 =[form_Energy_DF_total['Formation enegy (eV) (Na-rich)'][i] for i in range(31)]
array3 = [form_Energy_DF_total['Charge state'][i] for i in range(31)]
form_Energy_DF_combined=pd.DataFrame(index=index)
form_Energy_DF_combined['Charge state']=array3
form_Energy_DF_combined['Formation enegy (eV) (Na-rich)']=array2
form_Energy_DF_combined['Formation enegy (eV) (Na-poor)']=array1
#form_Energy_DF_combined[['Formation enegy (eV) (Na-poor)','Formation enegy (eV) (Na-rich)','Charge state']]=[array1,array2,array3]
#form_Energy_DF_combined
form_Energy_DF_total


# In[586]:
# Given the defect name, return charge state of the defect
# the key in dictionary Formation follows:
# V_Na^- or V_Mn^2-
# from this we can get the Charge state
def ObtChargeState(s):
    charge_str=s.split('^')[1]
    sign = 1
    if len(charge_str)==1:
        if charge_str == '0':
            charge_state = 0
        else:
            charge_state = 1
    else:
        charge_state = int(charge_str[0])
    if charge_str[-1] == '-':
        sign = -1
    if charge_str[-1] == '+':
        sign = 1
    return sign*charge_state


# In[715]:

# group same defects with different charge states
# return a dictionary: key=defects, values=(charge state, formation energy)
def group_Defect(chem_pot,elem):
    if elem.lower() == 'all':
        formation_energy = comp_energy_all(chem_pot)
    else:
        formation_energy = comp_energy_group(chem_pot,elem)
    defect_energy_by_group={}
    for i in formation_energy:
        defect_key = i.split('^')[0]
        defect_charge_state = ObtChargeState(i)
        if defect_key not in defect_energy_by_group:
            defect_energy_by_group[defect_key]=[(defect_charge_state,formation_energy[i])]
        else:
            defect_energy_by_group[defect_key].append((defect_charge_state,formation_energy[i]))
    for i in defect_energy_by_group:
        defect_energy_by_group[i].sort()
    return defect_energy_by_group


# In[719]:

#group_Defect(ChemPot1,'all')


# In[701]:

# calculate the charge state transition levels
# take keys of theoutput from group_Defect as input
# return a dictionary with transition levels
def transition_level_complicated(defect_string, grouped_defect_formation_energy):
    tran_lev_dict={}
    no_level = 'No transition level'
    if defect_string not in grouped_defect_formation_energy:
        return None
    number_of_charge_state = len(grouped_defect_formation_energy[defect_string])
    charge_state_list = grouped_defect_formation_energy[defect_string]
    key_string = ''
    if number_of_charge_state == 1:
        key_string = 'No transition level'
        tran_lev_dict[key_string] = 0
    for i in range(number_of_charge_state-1):
        key_string = '('+str(charge_state_list[i+1][0])+'/'+str(charge_state_list[i][0])+')'
        tran_lev_dict[key_string]= charge_state_list[i][1] - charge_state_list[i+1][1]
    return tran_lev_dict


# In[728]:

# calculate the charge state transition levels
# take keys of theoutput from group_Defect as input
# return a dictionary with transition levels
# no grouped_defect_formation_energy needed
# different from up
# group
def transition_level(defect_string):
    tran_lev_dict={}
    no_level = 'No transition level'
    global ChemPot1
    grouped_defect_formation_energy = group_Defect(ChemPot1,'all')
    if defect_string not in grouped_defect_formation_energy:
        return None
    number_of_charge_state = len(grouped_defect_formation_energy[defect_string])
    charge_state_list = grouped_defect_formation_energy[defect_string]
    key_string = ''
    if number_of_charge_state == 1:
        key_string = 'No transition level'
        tran_lev_dict[key_string] = 0
    for i in range(number_of_charge_state-1):
        key_string = '('+str(charge_state_list[i+1][0])+'/'+str(charge_state_list[i][0])+')'
        tran_lev_dict[key_string]= charge_state_list[i][1] - charge_state_list[i+1][1]
    return tran_lev_dict


####################################################################
# 
# set up transition kinks for plot
# (0, formation energy)
# (transition level, formation energy)
# (..., ...)
# (3.754, formation energy)
#
####################################################################
def set_up_plot(defect_type,chem_pot):
    #formation_energy=comp_energy_all(chem_pot)
    grouped_defect=group_Defect(chem_pot,'all')
    trans_points = [0]
    trans_points_energy = [grouped_defect[defect_type][-1][1]]
    trans_lev = transition_level(defect_type)
    if 'No transition level' not in trans_lev:
        for trans_key in trans_lev:
            trans_points.append(trans_lev[trans_key])
    trans_points.append(CBM)
    for i in range(1,len(trans_points)):
        trans_points_energy.append(grouped_defect[defect_type][-i][0]*trans_points[i]+grouped_defect[defect_type][-i][1])
    #trans_points.append(grouped_defect[defect_type][0][0]*CBM+grouped_defect[defect_type][0][1])
    # return the sorted list
    return sorted(list(zip(trans_points,trans_points_energy)))


# In[787]:

set_up_plot('Ni_Mn',ChemPot1)


####################################################################
# 
# print out the calculated transition levels:
# show the transition level of selected defect type
#
####################################################################
def print_tran_level(defect_group):
    temp = {}
    for key1 in defect_group:
        temp=transition_level(key1,defect_group)  
        print('The transition level of ', key1, ' is:\n')
        for key2 in temp:
            print(key2,': ', temp[key2],'\n')


####################################################################
# This function is used to calculate the binding energy of a polaron
# The inputs are: defect type 1, defect type 2, and polaron type 
# (whether it is a hole polaron or electron polaron)
# return a float-type value: binding energy (ev)
####################################################################
def comp_binding_polaron(defect1, defect2, polaron_type):
    # V_Na^0 and V_Na^- for example
    # polaron_type : hole
    global data_Energy, ChemPot1
    binding_energy = 0.0
    charge_state1 = ObtChargeState(defect1)
    charge_state2 = ObtChargeState(defect2)
    if abs(charge_state1-charge_state2)!=1:
        print('Please input two defect states with charge difference 1!')
        return None
    if polaron_type.lower() == 'hole':
        binding_energy = abs(comp_formation_energy(data_Energy,ChemPot1,defect1)-comp_formation_energy(data_Energy,ChemPot1,defect2))-hole_trap
    elif polaron_type.lower() == 'electron':
        binding_energy = CBM-abs(comp_formation_energy(data_Energy,ChemPot1,defect1)-comp_formation_energy(data_Energy,ChemPot1,defect2))-elct_trap
    return binding_energy
