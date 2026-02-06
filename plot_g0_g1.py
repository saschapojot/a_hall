import numpy as np
import matplotlib.pyplot as plt


#general case, a0 is not 0

a0_general=1
hbar=1
m=1

vR_list=[2,1,0.1,-0,1,-1,-2]

#edge case, a0=0
a0_edge=0

def g0(rho,a0,vR):
    return hbar**2/(2*m)*rho**2+(a0**2+vR**2*rho**2)**(1/2)

def g1(rho,a0,vR):
    return hbar ** 2 / (2 * m) * rho ** 2 - (a0 ** 2 + vR ** 2 * rho ** 2) ** (1 / 2)


# --- Plotting Logic ---

# Define the range for rho (momentum magnitude)
# Based on the document, rho >= 0. We'll plot from 0 to 4 to see the behavior clearly.
rho = np.linspace(0, 2, 400)
# Create a figure with two subplots: one for the General Case, one for the Edge Case
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# --- Plot 1: General Case (a0 = 1) ---
ax1.set_title(f"General Case ($a_0={a0_general}$)")
ax1.set_xlabel(r"$\rho$")
ax1.set_ylabel("Energy")
ax1.grid(True, linestyle='--', alpha=0.6)

for vR in vR_list:
    # Calculate values
    y_g0 = g0(rho, a0_general, vR)
    y_g1 = g1(rho, a0_general, vR)

    # Plot lines
    # We use the same color for g0 and g1 of the same vR, but different styles
    p = ax1.plot(rho, y_g0, label=f'$g_0, v_R={vR}$')
    color = p[0].get_color()
    ax1.plot(rho, y_g1, linestyle='--', color=color, label=f'$g_1, v_R={vR}$')

# Add legend to the first plot
ax1.legend()

# --- Plot 2: Edge Case (a0 = 0) ---
ax2.set_title(f"Edge Case ($a_0={a0_edge}$)")
ax2.set_xlabel(r"$\rho$")
ax2.set_ylabel("Energy")
ax2.grid(True, linestyle='--', alpha=0.6)

for vR in vR_list:
    # Calculate values
    y_g0 = g0(rho, a0_edge, vR)
    y_g1 = g1(rho, a0_edge, vR)

    # Plot lines
    p = ax2.plot(rho, y_g0, label=f'$g_0, v_R={vR}$')
    color = p[0].get_color()
    ax2.plot(rho, y_g1, linestyle='--', color=color, label=f'$g_1, v_R={vR}$')

# Add legend to the second plot
ax2.legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.savefig("g0g1.png")
plt.close()