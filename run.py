import openmc
import os
import subprocess

density_loadings = [0, 0.01, 0.1, 1, 5, 10, 25, 50, 75, 100, 150, 250, 500]

for loading in density_loadings: 
    name = f"{loading:06.2f}kgm3"
    run_dir = f"/Users/gli/code/FissileDependence/openmc_runs/{name}"
    print(f"Running OpenMC for {loading} kg/m3")
    
    subprocess.run(
        ["openmc"],
        cwd=run_dir,
        check=True
    )
