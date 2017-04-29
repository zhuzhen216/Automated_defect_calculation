
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
VBM = 0
CBM = 3.754
absolute_VBM = 3.014
# polaron self-trapping energy
elec_trap = 0.75
hole_trap = 0.86


# In[428]:

get_defect_info('V_Mn^')


# In[781]:

# set up a dictionary to defect_type, total_energy, correction
# 'Ni_Mn^+':[-1156.15,0.15,1]
data_Energy = {'bulk':[-1161.146,0.,0],'Na_Mn^0':[-1146.095,0.,0],'Mg_Na-V_Na^0':[-1159.163,0.0,0],'V_Na^-':[-1152.645,0.234,-1],'V_Na^0':[-1156.714,0.0,0],'Ni_Mn^-':[-1147.804,0.178,-1],'Ni_Mn^0':[-1152.156,0.,0],'Ge_Mn^+':[-1157.57,0.157,1],'Ge_Mn^0':[-1151.841,0.,0],'Ti_Mn^+':[-1165.538,0.134,1],'Ti_Mn^0':[-1159.923,0.,0],'Mg_Na^+':[-1167.622,0.16,1],'Mg_Na^0':[-1161.859,0.0,0],'Mg_Mn^-':[-1145.615,0.208,-1],'Mg_Mn^0':[-1149.804,0.,0],'Fe_Mn^0':[-1158.784,0.,0],'Co_Mn^0':[-1155.234,0.,0],'Al_Mn^0':[-1154.779,0.,0],'Mn_Na^0':[-1171.743,0.,0],'Mn_Na^+':[-1177.607,0.185,1],'Mn_Na^2+':[-1181.664,0.697,2],'Mn_Na-Na_Mn^0':[-1158.227,0.,0],'Mn_Na-Na_Mn^-':[-1153.676,0.187,-1],'V_Mn^3-':[-1127.904,1.879,-3],'V_Mn^2-':[-1132.523,0.984,-2],'V_Mn^-':[-1137.128,0.355,-1],'V_Mn^0':[-1140.706,0.,0],'H_i^0':[-1163.966,0.,0],'H_i^+':[-1169.478,0.22,1],'V_O^2+':[-1159.304,0.578,2],'V_O^+':[-1154.655,0.076,1],'V_O^0':[-1149.931,0.,0],'Na_Mn^2-':[-1137.182,0.776,-2],'Na_Mn^-':[-1141.656,0.265,-1]}


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

get_defect_info('V_Na-Mn_Na^0')


# In[489]:

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

group_Defect(ChemPot1,'all')


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


# In[786]:

transition_level('Ni_Mn')


# In[777]:

# set up transition kinks for plot
# (0, formation energy)
# (transition level, formation energy)
# ...
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


# In[721]:

# print out the calculated transition levels:
def print_tran_level(defect_group):
    temp = {}
    for key1 in defect_group:
        temp=transition_level(key1,defect_group)  
        print('The transition level of ', key1, ' is:\n')
        for key2 in temp:
            print(key2,': ', temp[key2],'\n')


# In[705]:

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


# In[479]:

comp_binding_polaron('V_Na^0','V_Na^-','hole')


# In[ ]:




# ---
# ## Al impurities
# ---

# In[368]:

# Al impurities
# at VBM
# under Na-rich conditions
Formation_Al_1 = {}
Formation_Al_1['Al_Na^2+']=-2.3715
Formation_Al_1['Al_Na-V_Na^+']=0.978
Formation_Al_1['Al_Na-V_2Na^0']=4.7475
Formation_Al_1['Al_Mn^0']=-0.3645
Formation_Al_1_CBM = {}
for i in Formation_Al_1:
    Formation_Al_1_CBM[i]=Formation_Al_1[i]+ObtChargeState(i)*CBM
Formation_Al_1_CBM


# In[514]:

comp_energy_group(ChemPot1,'Al')


# In[365]:

# at VBM
# under Na-rich conditions
Formation_Al_2 = {}
Formation_Al_2['Al_Na^2+']=-0.3085
Formation_Al_2['Al_Na-V_Na^+']=0.469
Formation_Al_2['Al_Na-V_2Na^0']=1.6665
Formation_Al_2['Al_Mn^0']=0.6635
Formation_Al_2_CBM = {}
for i in Formation_Al_1:
    Formation_Al_2_CBM[i]=Formation_Al_2[i]+ObtChargeState(i)*CBM
Formation_Al_2_CBM


# In[696]:

x = [VBM, CBM]
y1 = [[Formation_Al_1[i],Formation_Al_1_CBM[i]] for i in Formation_Al_1]
#print(y1)
#linetypes = ['solid','--',':','-.']
#fig = plt.figure()
for i in range(len(y1)):
    plt.plot(x,y1[i],linewidth=2.0)
plt.xlim((VBM,CBM))
plt.ylim((-4,6))
plt.xlabel('Fermi level (eV)')
plt.ylabel('Formation energy (eV)')
plt.xticks((0,1,2,3))
plt.axes().set_aspect(0.7)


# In[697]:

x = [VBM, CBM]
y2 = [[Formation_Al_2[i],Formation_Al_2_CBM[i]] for i in Formation_Al_2]
#linetypes = ['solid','--',':']
#fig = plt.figure()
for i in range(len(y1)):
    plt.plot(x,y2[i],linewidth=2.0)
plt.xlim((VBM,CBM))
plt.ylim((-4,6))
plt.xlabel('Fermi level (eV)')
plt.ylabel('Formation energy (eV)')
plt.xticks((0,1,2,3))
plt.axes().set_aspect(0.7)


# In[698]:

# try subplots for two sets of results with different chemical potential
x = [VBM, CBM]
y1 = [[Formation_Al_1[i],Formation_Al_1_CBM[i]] for i in Formation_Al_1]
y2 = [[Formation_Al_2[i],Formation_Al_2_CBM[i]] for i in Formation_Al_2]
f, ax = plt.subplots(1, 2, sharey=True)
for i in range(len(y1)):
    ax[0].plot(x, y1[i],linewidth=2)
    ax[1].plot(x, y2[i],linewidth=2)
for i in [0,1]:
    ax[i].set_xlim([0, 3.754])
    ax[i].set_xticks((0,1,2,3))
    ax[i].set_xlabel('Fermi level (eV)')
ax[0].set_ylabel('Formation energy (eV)')


# ---
# 1. Under Na-poor condition, Al_Na-V_2Na is the most stable defect for a large energy range within the gap
# 2. Al_Na can complex with V_Na to form defect complex
# 3. Under Na-rich condition, if the Fermi level is pinned near the middle of the gap, Al_Mn has lower formation energy than Al_Na, indicating that Al in Na layer can be effectively prohibited.
# ---

# ---
# ## Mg impurities
# ---

# In[499]:

# at VBM
# under Na-rich condition
#Formation_Mg_1={'Mg_Na^0':1.245,'Mg_Na^+':-1.344,'Mg_Na-V_Na^0':2.4435,'Mg_Mn^-':3.008,'Mg_Mn^0':1.625}
#Formation_Mg_2={'Mg_Na^0':1.763,'Mg_Na^+':-0.826,'Mg_Na-V_Na^0':0.3895,'Mg_Mn^-':2.491,'Mg_Mn^0':1.108}


# In[498]:

#defect_Mg=group_Defect(Formation_Mg_1)
#defect_Mg


# In[753]:

defect_Mg = group_Defect(ChemPot1, 'Mg')
defect_Mg


# In[754]:

Mg_Na_trans = transition_level('Mg_Na')
Mg_Na_trans


# In[196]:

#print_tran_level(defect_Mg)


# In[755]:

defect_Mg


# In[756]:

s_Na_rich = [set_up_plot(defect_key,ChemPot1) for defect_key in group_Defect(ChemPot1, 'Mg')]
x = [[s_Na_rich[i][j][0] for j in range(len(s_Na_rich[i]))] for i in range(len(s_Na_rich))]
y = [[s_Na_rich[i][j][1] for j in range(len(s_Na_rich[i]))] for i in range(len(s_Na_rich))]
linetypes = ['solid','--',':']
linecolor = ['r','g','b']
for i in range(len(s_Na_rich)):
    plt.plot(x[i],y[i],ls=linetypes[i],linewidth=2.0)
plt.ylim((-2,4))
plt.xlim((VBM,CBM))
plt.ylabel('Formation energy (eV)')
plt.xlabel('Fermi level (eV)')
plt.axes().set_aspect('equal')
plt.xticks(range(0,4))


# In[758]:

s_Na_poor = [set_up_plot(defect_key,ChemPot2) for defect_key in group_Defect(ChemPot2, 'Mg')]
x2 = [[s_Na_poor[i][j][0] for j in range(len(s_Na_poor[i]))] for i in range(len(s_Na_poor))]
y2 = [[s_Na_poor[i][j][1] for j in range(len(s_Na_poor[i]))] for i in range(len(s_Na_poor))]
linetypes = ['solid','--',':']
linecolor = ['r','g','b']
for i in range(len(s_Na_poor)):
    plt.plot(x2[i],y2[i],ls=linetypes[i],linewidth=2.0)
plt.ylim((-2,4))
plt.xlim((VBM,CBM))
plt.ylabel('Formation energy (eV)')
plt.xlabel('Fermi level (eV)')
plt.axes().set_aspect('equal')
#plt.xticks(,minor=True)
xtics = list(range(0,4))
#xtics_char = ['One','Two','Three','Four']
#ml = MultipleLocator(5)
plt.xticks(xtics)
#plt.axes().yaxis.set_minor_locator(ml)
#plt.xticks


# ---
# ## Mn-related defects
# ---

# In[312]:

#Formation_Mn_1={'Mn_Na^0':1.078,'Mn_Na^+':-1.587,'V_Mn^0':7.2675,'V_Mn^-':8.1865,'V_Mn^2-':10.4065,'V_Mn^3-':12.9065}


# In[793]:

Mn_plot_Na_rich = [set_up_plot(defect_key,ChemPot1) for defect_key in group_Defect(ChemPot1,'Mn')]
Mn_plot_Na_poor = [set_up_plot(defect_key,ChemPot2) for defect_key in group_Defect(ChemPot2,'Mn')]


# In[799]:

set_up_plot('V_Mn',ChemPot2)


# In[801]:

Mn_plot_Na_poor


# In[796]:

group_Defect(ChemPot1,'Mn')


# In[806]:

for list_i in Mn_plot_Na_poor:
    x_Mn_poor = [list_i[i][0] for i in range(len(list_i))]
    y_Mn_poor = [list_i[i][1] for i in range(len(list_i))]
    plt.plot(x_Mn_poor,y_Mn_poor,linewidth=2.0)
plt.xlim(VBM,CBM)
plt.ylim((-2,4))
plt.xticks([0,1,2,3])
plt.axes().set_aspect(1.2)
plt.xlabel('Fermi level (eV)')
plt.ylabel('Formation energy (ev)')


# In[802]:

for list_i in Mn_plot_Na_rich:
    x_Mn_rich = [list_i[i][0] for i in range(len(list_i))]
    y_Mn_rich = [list_i[i][1] for i in range(len(list_i))]
    plt.plot(x_Mn_rich,y_Mn_rich,linewidth=2.0)
plt.xlim(VBM,CBM)
plt.xticks([0,1,2,3])
plt.axes().set_aspect(0.6)
plt.xlabel('Fermi level (eV)')
plt.ylabel('Formation energy (ev)')


# In[ ]:



