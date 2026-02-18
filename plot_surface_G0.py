import numpy as np
import plotly.graph_objects as go

# --- 1. Physics Parameters ---
hbar = 1.0
m = 1.0
ma = 2.0
vR = 2.0
a0 = 0.5
theta_a = 0.1 * np.pi

# --- 2. Define the Function G0 ---
def calculate_G0(rho, gamma):
    # Term 1: Kinetic energy part ~ rho^2
    term1 = (hbar ** 2 / (2 * m)) * rho ** 2
    # Square root term
    sqrt_term = np.sqrt(a0 ** 2 + vR ** 2 * rho ** 2)
    # Term 2: Potential well
    term2 = sqrt_term

    # Term 3: Perturbation
    prefactor = (a0 * hbar ** 2) / (2 * ma)
    fraction = rho ** 2 / sqrt_term
    cosine_term = np.cos(2 * gamma - theta_a)
    term3 = prefactor * fraction * cosine_term
    return term1+term2+term3


resolution = 100
# Rho range
rho_range = np.linspace(0, 0.5, resolution)
# Gamma range: 0 to 2*pi
angle_range = np.linspace(0, 2 * np.pi, resolution)

# Create 2D Meshgrids
R, G = np.meshgrid(rho_range, angle_range)
Z = calculate_G0(R, G)

# Initialize Figure
fig = go.Figure()

# --- 3. Add Surface Trace ---
# Plotting directly on Rho (x-axis) and Gamma (y-axis)
fig.add_trace(go.Surface(
    z=Z,
    x=R,
    y=G,
    colorscale='Plasma',
    colorbar=dict(title='G0 Value')
))

# --- 4. Update Layout ---
fig.update_layout(
    title='G0 Potential Surface',
    scene=dict(
        xaxis_title='Rho (ρ)',
        yaxis_title='Gamma (γ)',
        zaxis_title='Energy G0(ρ, γ)',
        aspectmode='cube'
    ),
    width=900,
    height=800,
    margin=dict(l=65, r=50, b=65, t=90)
)

# --- 5. Save as Interactive HTML ---
output_filename = "G0_surface_plot.html"
fig.write_html(output_filename)

print(f"Plot saved successfully as '{output_filename}'")