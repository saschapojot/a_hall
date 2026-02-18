import numpy as np
import plotly.graph_objects as go

# --- 1. Physics Parameters (Common to both) ---
hbar = 1.0
m = 1.0
ma = 2.0
vR = 2.0
a0 = 0.5
theta_a = 0.1 * np.pi
resolution = 100

# --- 2. Define Functions ---

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
    return term1 + term2 + term3

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

# --- 3. Generate Data for G0 ---
# Range from script 1: rho 0 to 0.5
rho_range_0 = np.linspace(0, 0.5, resolution)
angle_range_0 = np.linspace(0, 2 * np.pi, resolution)
R0, G0 = np.meshgrid(rho_range_0, angle_range_0)
Z0 = calculate_G0(R0, G0)

# --- 4. Generate Data for G1 ---
# Calculate Critical Point Locations to determine range (from script 2)
n0, n1, n2 = 0, 1, 2
rho10 = np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) / (np.abs(vR) * hbar**2)

# First order corrections
rho11_n0 = a0 / (2 * vR**3) * (m**2 * vR**4 + a0**2 * hbar**4) / np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) * np.sign(vR) * (-1)**n0
rho11_n1 = a0 / (2 * vR**3) * (m**2 * vR**4 + a0**2 * hbar**4) / np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) * np.sign(vR) * (-1)**n1
rho11_n2 = a0 / (2 * vR**3) * (m**2 * vR**4 + a0**2 * hbar**4) / np.sqrt(m**2 * vR**4 - a0**2 * hbar**4) * np.sign(vR) * (-1)**n2

rho_c_n0 = rho10 + 1/ma * rho11_n0  # min
rho_c_n1 = rho10 + 1/ma * rho11_n1  # saddle
rho_c_n2 = rho10 + 1/ma * rho11_n2  # min

# Angular Locations
gamma_min_n0 = 0.5 * theta_a + 0.5 * n0 * np.pi
gamma_min_n2 = 0.5 * theta_a + 0.5 * n2 * np.pi
gamma_saddle_n1 = 0.5 * theta_a + 0.5 * n1 * np.pi

# Calculate Z values for markers
z_min_n0 = calculate_G1(rho_c_n0, gamma_min_n0)
z_min_n2 = calculate_G1(rho_c_n2, gamma_min_n2)
z_saddle_n1 = calculate_G1(rho_c_n1, gamma_saddle_n1)

# Range from script 2: rho 0 to 1.5 * rho_max
rho_max = max(rho_c_n0, rho_c_n1, rho_c_n2)
rho_range_1 = np.linspace(0, 1.5 * rho_max, resolution)
angle_range_1 = np.linspace(0, 2 * np.pi, resolution)
R1, G1 = np.meshgrid(rho_range_1, angle_range_1)
Z1 = calculate_G1(R1, G1)

# --- 5. Plotting ---
fig = go.Figure()

# Add Surface G0
fig.add_trace(go.Surface(
    z=Z0,
    x=R0,
    y=G0,
    colorscale='Plasma',
    opacity=0.85,
    name='G0',
    colorbar=dict(title='G0 Energy', x=1.02, len=0.7) # Position colorbar to the right
))

# Add Surface G1
fig.add_trace(go.Surface(
    z=Z1,
    x=R1,
    y=G1,
    colorscale='Viridis',
    opacity=0.85,
    name='G1',
    colorbar=dict(title='G1 Energy', x=1.15, len=0.7), # Position colorbar further right
    contours={
        "z": {"show": True, "usecolormap": True, "highlightcolor": "limegreen", "project": {"z": True}}
    }
))

# Add Markers for G1 Critical Points
# Marker: Minimum 1 (n=0)
fig.add_trace(go.Scatter3d(
    x=[rho_c_n0], y=[gamma_min_n0], z=[z_min_n0],
    mode='markers+text',
    marker=dict(size=5, color='red', symbol='circle'),
    text=["Min"], textposition="top center", textfont=dict(color='red'),
    name='G1 Min', showlegend=False
))

# Marker: Saddle Point (n=1)
fig.add_trace(go.Scatter3d(
    x=[rho_c_n1], y=[gamma_saddle_n1], z=[z_saddle_n1],
    mode='markers+text',
    marker=dict(size=5, color='orange', symbol='diamond'),
    text=["Saddle"], textposition="top center", textfont=dict(color='orange'),
    name='G1 Saddle', showlegend=False
))

# Marker: Minimum 2 (n=2)
fig.add_trace(go.Scatter3d(
    x=[rho_c_n2], y=[gamma_min_n2], z=[z_min_n2],
    mode='markers+text',
    marker=dict(size=5, color='red', symbol='circle'),
    text=["Min"], textposition="top center", textfont=dict(color='red'),
    name='G1 Min', showlegend=False
))

# --- 6. Layout Configuration ---
fig.update_layout(
    title='Combined Potential Surfaces: G0 and G1',
    width=1200,
    height=800,
    scene=dict(
        xaxis_title='Rho (ρ)',
        yaxis_title='Gamma (γ)',
        zaxis_title='Energy',
        yaxis=dict(
            range=[0, 2 * np.pi],
            tickvals=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
            ticktext=["0", "π/2", "π", "3π/2", "2π"]
        ),
        camera=dict(eye=dict(x=1.8, y=0.5, z=0.8)),
        aspectmode='cube'
    ),
    margin=dict(l=0, r=0, b=0, t=50)
)

# --- 7. Save as Interactive HTML ---
output_filename = "Combined_G0_G1_Surface.html"
fig.write_html(output_filename)
print(f"Combined plot saved successfully to '{output_filename}'")