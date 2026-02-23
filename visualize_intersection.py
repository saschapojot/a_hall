import numpy as np
import plotly.graph_objects as go

def visualize_intersection_interactive():
    # 1. Configuration
    a = 4.0               # The value for t0^2 + t1^2 = a
    radius = np.sqrt(a)   # Radius of the intersection circle

    # 2. Create the Surface (Paraboloid)
    limit = radius * 1.5
    grid_size = 100
    t0 = np.linspace(-limit, limit, grid_size)
    t1 = np.linspace(-limit, limit, grid_size)
    T0, T1 = np.meshgrid(t0, t1)
    Z = T0**2 + T1**2

    # 3. Create the Intersection Curve (Circle at height z=a)
    theta = np.linspace(0, 2 * np.pi, 200)
    t0_curve = radius * np.cos(theta)
    t1_curve = radius * np.sin(theta)
    z_curve = np.full_like(t0_curve, a)

    # 4. Calculate Gradient Vectors (Cones)
    num_vectors = 16
    theta_vec = np.linspace(0, 2 * np.pi, num_vectors, endpoint=False)
    t0_vec = radius * np.cos(theta_vec)
    t1_vec = radius * np.sin(theta_vec)
    z_vec = np.full_like(t0_vec, a)

    # Gradient of h = t0^2 + t1^2 is <2t0, 2t1>
    u = 2 * t0_vec
    v = 2 * t1_vec
    w = np.zeros_like(u) # Keep vectors flat on the z-plane

    # Normalize vectors for consistent cone size
    norm = np.sqrt(u**2 + v**2 + w**2)
    u, v, w = (u/norm, v/norm, w/norm)

    # 5. Build the Plotly Figure
    fig = go.Figure()

    # Add Surface
    fig.add_trace(go.Surface(
        z=Z, x=T0, y=T1,
        colorscale='Viridis',
        opacity=0.6,
        name='Surface z=t₀²+t₁²'
    ))

    # Add Intersection Curve
    fig.add_trace(go.Scatter3d(
        x=t0_curve, y=t1_curve, z=z_curve,
        mode='lines',
        line=dict(color='red', width=5),
        name=f'Intersection (z={a})'
    ))

    # Add Gradient Vectors
    fig.add_trace(go.Cone(
        x=t0_vec, y=t1_vec, z=z_vec,
        u=u, v=v, w=w,
        sizemode="absolute",
        sizeref=0.5,
        anchor="tail",
        colorscale=[[0, 'black'], [1, 'black']],
        showscale=False,
        name='Gradient Direction'
    ))

    # 6. Layout Settings with renamed axes
    fig.update_layout(
        title=f'Interactive Intersection: z = t₀² + t₁² at z={a}',
        scene=dict(
            xaxis_title='t₀',
            yaxis_title='t₁',
            zaxis_title='Z',
            aspectmode='cube'
        ),
        margin=dict(l=0, r=0, b=0, t=50)
    )

    # 7. Save to HTML
    output_file = "intersection_t0_t1.html"
    fig.write_html(output_file)
    print(f"Interactive plot saved successfully to {output_file}")

if __name__ == "__main__":
    try:
        visualize_intersection_interactive()
    except ImportError:
        print("Error: This script requires 'plotly'.")
        print("Please install it using: pip install plotly numpy")