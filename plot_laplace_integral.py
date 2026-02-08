import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return -x**2+1

x_range_0=np.linspace(0,1,100)
x_range_1=np.linspace(0.5,1,100)

y_range_0=f(x_range_0)
y_range_1=f(x_range_1)
# Create a figure with two subplots side by side (1 row, 2 columns)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
# Plot 1: Full Range (0 to 1)
axes[0].plot(x_range_0, y_range_0, color='blue', label='Range [0, 1]')
axes[0].set_ylabel('h(x)')
axes[0].set_title(r"$h^{'}(a)=0$")
# Remove x-ticks and label x=0 as 'a'
axes[0].set_xticks([])  # Remove standard ticks
axes[0].set_yticks([])  # Remove standard ticks
axes[0].text(0, -0.15, 'a', ha='center', va='top', fontsize=12, color='black') # Place 'a' below 0
axes[0].axvline(x=0, color='black', linestyle='--', linewidth=1)
# Add Figure Label (a)
# transform=axes[0].transAxes uses coordinates relative to the box (0,0 is bottom-left, 1,1 is top-right)
axes[0].text(0.5, -0.2, '(a)', transform=axes[0].transAxes,
             ha='center', va='top', fontsize=14, color='black')
# Label 'b' at x=1 (last value)
axes[0].text(1, -0.15, 'b', ha='center', va='top', fontsize=12, color='black')
axes[0].axvline(x=1, color='black', linestyle='--', linewidth=1)

# Plot 2: Partial Range (0.5 to 1)
axes[1].plot(x_range_1, y_range_1, color='red', label='Range [0.5, 1]')
axes[1].set_xlabel('x')
axes[1].set_ylabel('h(x)')
axes[1].grid(True)
axes[1].set_xticks([]) # Remove standard ticks
axes[1].set_yticks([])  # Remove standard ticks
axes[1].text(0.5, -0.15, 'a', ha='center', va='top', fontsize=12, color='black')
axes[1].axvline(x=0.5, color='black', linestyle='--', linewidth=1)
axes[1].set_title(r"$h^{'}(a)<0$")
# Label 'b' at x=1 (last value)
axes[1].text(1, -0.15, 'b', ha='center', va='top', fontsize=12, color='black')
axes[1].axvline(x=1, color='black', linestyle='--', linewidth=1)
# Add Figure Label (b)
# Placed slightly lower (-0.25) to account for the 'x' label on this plot
axes[1].text(0.5, -0.25, '(b)', transform=axes[1].transAxes,
             ha='center', va='top', fontsize=14, color='black')
# Optional: Set consistent y-axis limits to make visual comparison easier
# The max y is 1 (at x=0) and min y is 0 (at x=1)
axes[0].set_ylim(-0.1, 1.1)
axes[1].set_ylim(-0.1, 1.1)
plt.tight_layout()
plt.savefig("laplace.png")
plt.close()



