import numpy as np
import plotly.graph_objects as go

# --- 1. Physics Parameters ---
hbar = 1.0
m = 1.0
ma = 2.0
vR = 2.0
a0 = 0.5
theta_a = 0.1 * np.pi


# --- 2. Define the Function G1 ---
def calculate_G1(rho, gamma):
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

    return term1 - term2 - term3


# --- 3. Calculate Critical Point Locations ---

# Radial Center (rho_c)
n0 = 0
n1 = 1
n2 = 2
rho10 = np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) / (np.abs(vR) * hbar**2)

# First order corrections
rho11_n0 = a0 / (2 * vR**3) * (m**2 * vR**4 + a0**2 * hbar**4) / np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) * np.sign(vR) * (-1)**n0
rho11_n1 = a0 / (2 * vR**3) * (m**2 * vR**4 + a0**2 * hbar**4) / np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) * np.sign(vR) * (-1)**n1
rho11_n2 = a0 / (2 * vR**3) * (m**2 * vR**4 + a0**2 * hbar**4) / np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) * np.sign(vR) * (-1)**n2

rho_c_n0 = rho10 + 1/ma * rho11_n0  # min
rho_c_n1 = rho10 + 1/ma * rho11_n1  # saddle
rho_c_n2 = rho10 + 1/ma * rho11_n2  # min

# Angular Locations based on user input:
gamma_min_n0 = 0.5 * theta_a + 0.5 * n0 * np.pi  # n=0
gamma_min_n2 = 0.5 * theta_a + 0.5 * n2 * np.pi  # n=2
gamma_saddle_n1 = 0.5 * theta_a + 0.5 * n1 * np.pi  # n=1

# Calculate Z values for markers
z_min_n0 = calculate_G1(rho_c_n0, gamma_min_n0)
z_min_n2 = calculate_G1(rho_c_n2, gamma_min_n2)
z_saddle_n1 = calculate_G1(rho_c_n1, gamma_saddle_n1)

# --- 4. Generate Mesh Data ---
resolution = 100
rho_max = max(rho_c_n0, rho_c_n1, rho_c_n2)

# Rho range
rho_range = np.linspace(0, 1.5 * rho_max, resolution)

# Gamma range: 0 to 2*pi
angle_range = np.linspace(0, 2 * np.pi, resolution)

# Create 2D Meshgrids
R, G = np.meshgrid(rho_range, angle_range)
Z = calculate_G1(R, G)

# Initialize Figure
fig = go.Figure()

# 1. Add the Surface Plot
fig.add_trace(go.Surface(
    z=Z,
    x=R,
    y=G,
    colorscale='Viridis',
    opacity=0.9,
    colorbar=dict(title='Energy'),
    contours={
        "z": {"show": True, "usecolormap": True, "highlightcolor": "limegreen", "project": {"z": True}}
    }
))

# 2. Add Marker: Minimum 1 (n=0)
fig.add_trace(go.Scatter3d(
    x=[rho_c_n0],
    y=[gamma_min_n0],
    z=[z_min_n0],
    mode='markers+text',
    marker=dict(size=6, color='red', symbol='circle'),
    text=["Min"],
    textposition="top center",
    textfont=dict(color='red'),  # MODIFICATION: Set text color to red
    name='Min (n=0)',
    showlegend=False
))

# 3. Add Marker: Saddle Point (n=1)
fig.add_trace(go.Scatter3d(
    x=[rho_c_n1],
    y=[gamma_saddle_n1],
    z=[z_saddle_n1],
    mode='markers+text',
    marker=dict(size=6, color='orange', symbol='diamond'),
    text=["Saddle"],
    textposition="top center",
    textfont=dict(color='orange'),  # MODIFICATION: Set text color to orange
    name='Saddle (n=1)',
    showlegend=False
))

# 4. Add Marker: Minimum 2 (n=2)
fig.add_trace(go.Scatter3d(
    x=[rho_c_n2],
    y=[gamma_min_n2],
    z=[z_min_n2],
    mode='markers+text',
    marker=dict(size=6, color='red', symbol='circle'),
    text=["Min"],
    textposition="top center",
    textfont=dict(color='red'),  # MODIFICATION: Set text color to red
    name='Min (n=2)',
    showlegend=False
))

# 5. Layout Configuration
fig.update_layout(
    title='G1 Potential Surface',
    width=1200,
    height=800,
    scene=dict(
        xaxis_title='Rho (ρ)',
        yaxis_title='Gamma (γ)',
        zaxis_title='Energy G1(ρ, γ)',
        yaxis=dict(
            range=[0, 2 * np.pi],
            tickvals=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
            ticktext=["0", "π/2", "π", "3π/2", "2π"]
        ),
        camera=dict(
            eye=dict(x=1.8, y=0.5, z=0.8)
        ),
        aspectratio=dict(x=1, y=1, z=0.6)
    ),
    margin=dict(l=0, r=0, b=0, t=50)
)

# Save and Show
output_filename = "G1_surface_topology.html"
fig.write_html(output_filename)
print(f"Plot saved to {output_filename}")