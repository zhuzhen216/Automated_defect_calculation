# defect_analysis
### about functions to do data analysis
#
### add_TolEnergy_correction_state(str defect_name, list [tot energy, correction, charge state])
### return type: None; just add to data_Energy dictionary
#
### get_defect_info(str defect_name)
### return type: dictionary {elem1:num, elem2:num,...}. This will relate to the formation energy calculation.
#
### comp_formation_energy(dict chem_pot, str defect_name)
### return a float number: the formation energy of the defect.
#
### comp_energy_group(dict chem_pot, str element_char)
### return a dictionary {str defect_name: float formation energy}
### to calculate the formation energy of a group of defects containing the same element
#
### comp_energy_all(dict chem_pot)
### return a dictionary {str defect_name: float formation energy}
### for all the defects
#
### ObtChargeState(str defect_name)
### return the charge state of the defect
