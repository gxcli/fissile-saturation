import numpy as np
import openmc
import matplotlib.pyplot as plt

density_loadings = [0.1, 5, 10, 50, 100, 500] # [0, 0.01, 0.1, 1, 5, 10, 25, 50, 75, 100, 150, 250, 500]
energy_bins = np.logspace(-1, 7.3, 100)
results = []

for loading in density_loadings: 
    name = f"{loading:06.2f}kgm3"
    run_dir = f"/Users/gli/code/FissileDependence/openmc_runs/{name}"
    print(f"Collecting tallies for {loading} kg/m3")
    
    sp = openmc.StatePoint(run_dir + "/statepoint.150.h5")
    tally = sp.get_tally(name="U238(n,gamma) to U239")
    result = tally.mean.flatten()
    results.append(result)
    
    cumulative_result = np.cumsum(result) / np.sum(result)
    plt.stairs(cumulative_result, energy_bins, label=f"{loading} kg/m3", linewidth=2)

results = np.array(results)
print(results.shape)

# PLOTTING ###################################################
# # cdf plot
plt.legend()
plt.xscale('log')
plt.xlabel('Neutron Energy (eV)')   
plt.ylabel(r'CDF of U238$(n,\gamma)$ to Pu239 Tallies')
plt.ylim(0,1)
plt.savefig('U238_ngamma_cdf.png', dpi=300)


# plt.figure()

# for i in range(10, 99, 10):
#     lower, upper = energy_bins[i], energy_bins[i + 1] 
#     print(lower, upper)
#     plt.scatter(density_loadings, results[:, i]/np.max(results[:, i]), label=f"{lower:.2e} - {upper:.2e} eV", s=20)
    
# plt.plot(np.linspace(0, 500, 100), np.linspace(0, 1, 100), linewidth=1, color='black', linestyle='--')
# plt.xlabel('Density Loading (kg/m3)', fontsize=12)
# plt.ylabel('U238(n,gamma) Tallies \n Normalized Against Max. Tally per Energy Bin', fontsize=12)
# plt.legend(fontsize=10)
# plt.savefig('U238_ngamma_vs_loading.png', dpi=300)