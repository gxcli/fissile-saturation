import numpy as np
import openmc 
import os

# 
# DOUBLE CHECK CROSS SECTIONS

# PARAMETERS #########################################################################################################
# natural lithium enrichment, uranium loaded, 900.0 K temperature 

u238_density_loading_kgm3 = 10 # kg/m3, of total breeder volume 
print(u238_density_loading_kgm3, "kg/m3 of U238 in breeder")
u238_density_loading_gcm3 = u238_density_loading_kgm3 / 1000 # g/cm3

name = f"{u238_density_loading_kgm3:06.2f}kgm3"
run_dir = f"./openmc_runs/{name}"
os.makedirs(run_dir, exist_ok=True)


# MATERIALS #########################################################################################################
# reference values: dt-fusion-illicit/Python/parameters.py
print('MATERIAL TIME!')

flibe = openmc.Material(name='2(LiF)-BeF2', temperature=900.0)
flibe.set_density('g/cm3', 1.9505) 
flibe.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='ao', enrichment=7.50)
# natural enrichment of lithium 

fertile = openmc.Material(name='fertile', temperature=900.0)
fertile.add_elements_from_formula('UF4', enrichment_type='ao')
uf4_density = 6.88 # g/cm3
fertile.set_density('g/cm3', uf4_density)

# UF4 displaces its own volume in FLiBe 
AMU_F19 = 18.9984   # g/mol 
AMU_U235 = 235.0439299  # g/mol 
AMU_U238 = 238.05078826 # g/mol 

U235_WTPERC_ENRICH = 0.71     # weight percentage
U238_WTPERC_ENRICH = 100 - U235_WTPERC_ENRICH     # weight percentage

# assuming 100g of uranium mixture with weight percentages above, 
u235_mol, u238_mol = U235_WTPERC_ENRICH / AMU_U235, U238_WTPERC_ENRICH / AMU_U238  # mols
total_mol = u235_mol + u238_mol
U235_MOLFRAC_ENRICH, U238_MOLFRAC_ENRICH = u235_mol / total_mol, u238_mol / total_mol # molar fractions

MASSFRAC_U238_OF_UF4 = (AMU_U238 * U238_MOLFRAC_ENRICH) / (AMU_U235 * U235_MOLFRAC_ENRICH + AMU_U238 * U238_MOLFRAC_ENRICH + 4 * AMU_F19)

u238_density_of_uf4 = (MASSFRAC_U238_OF_UF4 * uf4_density) # g U238 / cm3 
vf_uf4 = u238_density_loading_gcm3 / u238_density_of_uf4 # volume fraction of UF4 in breeder
vf_flibe = 1 - vf_uf4 # volume fraction of FLiBe in breeder
print("volume fraction of UF4 in FLiBe:", vf_uf4)

blanket = openmc.Material(name=f"{u238_density_loading_kgm3:07.2f} kg/m3 | {(vf_uf4*100):.4f} vol% of breeder'", temperature=900.0)
blanket = openmc.Material.mix_materials([flibe, fertile], [vf_flibe, vf_uf4], 'vo')

materials = openmc.Materials([blanket])
materials.export_to_xml(path=run_dir + '/materials.xml')



# check that the cross sections are correct 



# GEOMETRY #########################################################################################################
print('GEOMETRY TIME!')