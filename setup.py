import numpy as np
import openmc 
import os



# PARAMETERS #########################################################################################################
# natural lithium enrichment, uranium loaded, 900.0 K temperature
# cross sections - need to download ENDF VIII.0 version and set the global environment variable 

density_loadings = [0, 0.01, 0.1, 1, 5, 10, 25, 50, 75, 100, 150, 250, 500]

u238_density_loading_kgm3 =  500 # kg/m3, of total breeder volume 
print(u238_density_loading_kgm3, "kg/m3 of U238 in breeder")
u238_density_loading_gcm3 = u238_density_loading_kgm3 / 1000 # g/cm3

name = f"{u238_density_loading_kgm3:06.2f}kgm3"
run_dir = f"/Users/gli/code/FissileDependence/openmc_runs/{name}"
os.makedirs(run_dir, exist_ok=True)


# MATERIALS #########################################################################################################
# reference values: dt-fusion-illicit/Python/parameters.py
print('MATERIAL TIME!')

# breeder material
flibe = openmc.Material(name='2(LiF)-BeF2', temperature=900.0)
flibe.set_density('g/cm3', 1.9505) 
flibe.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='ao', enrichment=7.50)

# constants for calculating UF4
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

# UF4 displaces its own volume in FLiBe 
fertile = openmc.Material(name='fertile', temperature=900.0)
fertile.add_elements_from_formula('UF4', enrichment_type='ao') # atomic ratio

uf4_density = 6.88 # g/cm3
fertile.set_density('g/cm3', uf4_density) 

u238_density_of_uf4 = (MASSFRAC_U238_OF_UF4 * uf4_density) # g U238 / cm3 
vf_uf4 = u238_density_loading_gcm3 / u238_density_of_uf4 # volume fraction of UF4 in breeder
vf_flibe = 1 - vf_uf4 # volume fraction of FLiBe in breeder
print("volume fraction of UF4 in FLiBe:", vf_uf4)

blanket = openmc.Material.mix_materials([flibe, fertile], [vf_flibe, vf_uf4], 'vo')
blanket.name = f"{u238_density_loading_kgm3:07.2f} kg/m3 | {(vf_uf4*100):.4f} vol% of breeder"
blanket.temperature = 900.0

materials = openmc.Materials([blanket])
materials.export_to_xml(path=run_dir + '/materials.xml')

# GEOMETRY #########################################################################################################
print('GEOMETRY TIME!')
# 800x800x800 cm box, ballpark of the 773.3 m^3 blanket
print('volume of box = 512 m3')

# surfaces 
left, right = openmc.XPlane(x0=-400, name='left', boundary_type='periodic'), openmc.XPlane(x0=400, name='right', boundary_type='periodic')
bottom, top = openmc.YPlane(y0=-400, name='bottom', boundary_type='periodic'), openmc.YPlane(y0=400, name='top', boundary_type='periodic')
back, front = openmc.ZPlane(z0=-400, name='back', boundary_type='periodic'), openmc.ZPlane(z0=400, name='front', boundary_type='periodic')

# boundary condition
back.periodic_surface = front 
bottom.periodic_surface = top
left.periodic_surface = right

# region 
box = +left & -right & + bottom & -top & +back & -front

# cells
cell = openmc.Cell(name='breeder', fill=blanket, region=box)

geom = openmc.Geometry([cell])
geom.export_to_xml(path=run_dir + '/geometry.xml')

# SETTINGS #########################################################################################################
print('SETTINGS TIME!')
settings = openmc.Settings()
settings.run_mode = 'fixed source' 

# source
source=openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.0e6], [1.0]) # 14 MeV source
settings.source = source

# run strat
settings.batches = 150
settings.particles = int(5e5)

# cut off 
settings.cutoff = {'energy_neutron': 1e-1}

# temperature check 
settings.temperature_default = 900.0

# output files 
settings.output = {'path': run_dir}
settings.export_to_xml(path=run_dir + '/settings.xml')


# TALLIES ##########################################################################################################
print('TALLIES TIME!')
tallies = openmc.Tallies()

# U238->Pu239 filters and scores
tally = openmc.Tally(name='U238(n,gamma) to Pu239')
energy_bins = np.logspace(-1, 7.3, 100) # energy bins from 0.1 eV to ~20 MeV
energy_filter = openmc.EnergyFilter(energy_bins) # energy bins or can this be done continuously? 
tally.filters = [energy_filter]
tally.scores = ['(n,gamma)'] # to produce Pu239
tally.nuclides = ['U238']

tallies.append(tally)
tallies.export_to_xml(path=run_dir + '/tallies.xml')
