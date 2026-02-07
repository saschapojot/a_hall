import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

#general case, a0 is not 0
#large |vR|
a0_general=1
hbar=1
m=1


#edge case, a0=0
a0_edge=0

def g0(rho,a0,vR):
    return hbar**2/(2*m)*rho**2+(a0**2+vR**2*rho**2)**(1/2)

def g1(rho,a0,vR):
    return hbar ** 2 / (2 * m) * rho ** 2 - (a0 ** 2 + vR ** 2 * rho ** 2) ** (1 / 2)


def rho_positive(vR,a0):
    val=np.sqrt(m**2*vR**4-a0**2*hbar**4)/(np.abs(vR)*hbar**2)
    return val


def plot_large_abs_vR(vR,a0):
    """

    :param vR: large |vR|
    :param a0: a0 not 0
    :return:
    """
    if not vR**2>np.abs(a0)*hbar**2/m:
        raise ValueError("vR**2>np.abs(a0)*hbar**2/m not satisfied")
    if np.abs(a0)==0:
        raise ValueError("a0 is 0")


    rho_solution=rho_positive(vR,a0)
    rho_range=np.linspace(0,2.5*rho_solution,300)
    out_dir="./large_abs_vR/"
    Path(out_dir).mkdir(exist_ok=True,parents=True)
    fig,ax=plt.subplots(1, 1, figsize=(5, 5))
    ax.set_title(r"General Case ($|a_0|\neq 0$), $v_{R}^{2}>\frac{|a_0|\hbar^2}{m}$")
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel("Energy")
    ax.grid(True, linestyle='--', alpha=0.6)
    # Remove Y-axis numerical labels
    ax.set_yticklabels([])
    y_g0 = g0(rho_range, a0, vR)
    y_g1 = g1(rho_range, a0, vR)
    ax.plot(rho_range,y_g0,color="blue",label="$g_0$")
    ax.plot(rho_range,y_g1,color="green",label="$g_1$")
    #add horizontal lines
    g1_min=-(m**2*vR**4+a0**2*hbar**4)/(2*m*vR**2*hbar**2)
    minus_abs_a0=-np.abs(a0)
    abs_a0=np.abs(a0)
    # Dictionary mapping values to their LaTeX labels
    special_vals = {
        abs_a0: r"$|a_0|$",
        minus_abs_a0: r"$-|a_0|$",
        g1_min: r"$(g_{1})_{\text{min}}$"
    }
    # Get the max x-value to position text on the right side

    for val, label in special_vals.items():
        # Draw the line
        ax.axhline(y=val, color='black', linestyle='--', linewidth=1)

        # Add the text label
        # x position is slightly past the plot edge (x_max * 1.02)
        # va='center' aligns the text vertically with the line
        ax.text(-1, val, f"  {label}", verticalalignment='center', fontsize=10)



    ax.legend()
    plt.tight_layout()
    plt.savefig(out_dir+f"/vR_{vR}_a0_{a0}.png")
    plt.savefig(out_dir + f"/vR_{vR}_a0_{a0}.svg")
    plt.close()




def plot_plot_large_abs_vR_with_muF(vR,a0):
    """

    :param vR: large |vR|
    :param a0: a0 not 0
    :return:
    """
    if not vR ** 2 > np.abs(a0) * hbar ** 2 / m:
        raise ValueError("vR**2>np.abs(a0)*hbar**2/m not satisfied")
    if np.abs(a0) == 0:
        raise ValueError("a0 is 0")
    rho_solution = rho_positive(vR, a0)
    rho_range = np.linspace(0, 2.7 * rho_solution, 300)
    out_dir = "./large_abs_vR/"
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    #  horizontal line values
    g1_min = -(m ** 2 * vR ** 4 + a0 ** 2 * hbar ** 4) / (2 * m * vR ** 2 * hbar ** 2)
    minus_abs_a0 = -np.abs(a0)
    abs_a0 = np.abs(a0)
    # Dictionary mapping values to their LaTeX labels
    special_vals = {
        abs_a0: r"$|a_0|$",
        minus_abs_a0: r"$-|a_0|$",
        g1_min: r"$(g_{1})_{\text{min}}$"
    }
    muF_list = [g1_min, (g1_min + minus_abs_a0) / 2, minus_abs_a0, (minus_abs_a0 + abs_a0) / 2, abs_a0, abs_a0 + 1]
    for counter, muF in enumerate(muF_list):
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.set_title(r"General Case ($|a_0|\neq 0$), $v_{R}^{2}>\frac{|a_0|\hbar^2}{m}$, " + fr"$\mu_{{F{counter}}}$")
        ax.set_xlabel(r"$\rho$")
        ax.set_ylabel("Energy")
        ax.grid(True, linestyle='--', alpha=0.6)
        # Remove Y-axis numerical labels
        ax.set_yticklabels([])
        y_g0 = g0(rho_range, a0, vR)
        y_g1 = g1(rho_range, a0, vR)
        ax.plot(rho_range, y_g0, color="blue", label="$g_0$")
        ax.plot(rho_range, y_g1, color="green", label="$g_1$")
        # add horizontal lines
        for val, label in special_vals.items():
            # Draw the line
            ax.axhline(y=val, color='black', linestyle='--', linewidth=0.7)

            # Add the text label
            # x position is slightly past the plot edge (x_max * 1.02)
            # va='center' aligns the text vertically with the line
            ax.text(-1, val, f"  {label}", verticalalignment='center', fontsize=10)
        ax.axhline(y=muF, color='red', linestyle='-', linewidth=1,label="chemical potential")
        ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir + f"/vR_{vR}_a0_{a0}_muF{counter}.png")
        plt.savefig(out_dir + f"/vR_{vR}_a0_{a0}_muF{counter}.svg")
        plt.close()



vR_large_list=[2]
for vR in vR_large_list:
    plot_large_abs_vR(vR,a0_general)
    plot_plot_large_abs_vR_with_muF(vR,a0_general)