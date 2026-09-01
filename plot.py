import numpy as np
import openmc
import matplotlib.pyplot as plt

# make the histogram
# plot just one of them 
# and then overlay them on top of each other

density_loadings = [0.1, 5, 10, 50, 100, 500] # [0, 0.01, 0.1, 1, 5, 10, 25, 50, 75, 100, 150, 250, 500]

# set up the plot bare necessities here
energy_bins = np.logspace(-1, 7.3, 100)


for loading in density_loadings: 
    name = f"{loading:06.2f}kgm3"
    run_dir = f"/Users/gli/code/FissileDependence/openmc_runs/{name}"
    print(f"Collecting tallies for {loading} kg/m3")
    
    sp = openmc.StatePoint(run_dir + "/statepoint.150.h5")
    tally = sp.get_tally(name="U238(n,gamma) to U239")
    results = tally.mean.flatten()
    
    # plt.stairs(results, energy_bins, label=f"{loading} kg/m3")
    cumulative_results = np.cumsum(results) / np.sum(results)
    plt.stairs(cumulative_results, energy_bins, label=f"{loading} kg/m3", linewidth=2)
    

plt.legend()
plt.xscale('log')
plt.xlabel('Neutron Energy (eV)')   
plt.ylabel('CDF of U238(n,gamma) to U239 Rxn')
plt.ylim(0,1)
plt.savefig('U238_ngamma_cdf.png', dpi=300)