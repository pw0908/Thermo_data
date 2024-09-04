import matplotlib.pyplot as plt
import numpy as np
import CoolProp
import matplotlib
import sys
import json

def graph_density_error(subs, k, CES):
    handle = CoolProp.AbstractState("HEOS", subs)
    pc = handle.p_critical()
    Tc = handle.T_critical()
    
    # Load metadata and data files once
    with open(f'{subs} metadata.json') as f:
        metadata = json.load(f)

    predicted_density_data = np.load(f'{CES} {subs} density.npz')
    experimental_density_data = np.load(f'Ld CoolProp {subs}.npz')
    VP_error_data = np.load(f'{CES}_{subs}_density_sat.npz')

    # Extract data arrays
    predicted_density = predicted_density_data["Phi_PR"]
    experimental_density = experimental_density_data["Phi_PR"]

    # Since predicted_density and experimental_density have different axes order,
    # we need to transpose predicted_density to match the axes of experimental_density
    predicted_density = np.transpose(predicted_density)

    Tvp = VP_error_data["Temperature"]
    VP = VP_error_data["VaporPressure"]
    error_lig = VP_error_data["Error_liq_sat"]
    errors_gas = VP_error_data["Error_gas_sat"]
    errors = (error_lig + errors_gas) / 2

    # Calculate the error
    Error = np.transpose(np.abs(predicted_density - experimental_density) * 100 / experimental_density)

    # Precompute meshgrid
    P, T = np.meshgrid(metadata["grid"]["P"], metadata["grid"]["T"])

    # Define common parameters for the plots
    levels = np.linspace(0, 30, 11)
    cmap = plt.get_cmap('RdYlGn_r')
    cmap.set_over('red')
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False)

    # Plot setup function to avoid code repetition
    def setup_plot():
        plt.contourf(T / Tc, P / pc, Error, levels=levels, cmap=cmap, extend='max')
        cbar = plt.colorbar(label='Error (%)')
        cbar.set_ticks(levels)
        cbar.set_ticklabels([str(int(lvl)) for lvl in levels])
        plt.grid()
        plt.axvline(x=1, linestyle='--', linewidth=3, color="k")
        plt.axhline(y=1, linestyle='--', linewidth=3, color="k")
        plt.plot(Tvp / Tc, VP / pc, linestyle='-', linewidth=3, color='k')
        scatter_colors = cmap(norm(errors))
        plt.scatter(Tvp / Tc, VP / pc, c=scatter_colors, edgecolor='black', s=50, zorder=2)

    # First figure
    plt.figure(2 ** k)
    plt.title(f"{CES} Density error for {subs}")
    plt.yscale("log")
    setup_plot()
    plt.xlabel("Tr")
    plt.ylabel("Pr")
    plt.gca().set_ylim(bottom=0.001)
    plt.savefig(f'{CES} {subs} density error big.png')
    plt.close()

    # Second figure
    plt.figure(3 ** k)
    plt.title(f"{CES} Density error for {subs}")
    setup_plot()
    plt.xlabel("Tr")
    plt.yscale("log")
    plt.ylabel("Pr")
    plt.xlim(0.5, 1.5)
    plt.ylim(0.001, 10)
    plt.savefig(f'{CES} {subs} density error small.png')
    plt.close()

# Call the function
graph_density_error(sys.argv[1], 3, sys.argv[2])
