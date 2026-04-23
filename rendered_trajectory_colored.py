from ovito.io import import_file
from ovito.modifiers import PythonScriptModifier
from ovito.vis import Viewport, OpenGLRenderer

# -----------------------------
# User settings
# -----------------------------
input_file = "trajectory.xyz"
output_file = "argon_movie.mp4"

# Rahman-style box length
L = 10.229343992174416

# -----------------------------
# Load trajectory
# -----------------------------
pipeline = import_file(input_file, multiple_frames=True)

# -----------------------------
# Add color/radius properties
# -----------------------------
def setup_particle_visuals(frame, data):
    import numpy as np

    positions = data.particles.positions
    n = len(positions)

    # Create per-particle color property if absent
    if "Color" not in data.particles:
        colors = data.particles_.create_property("Color", data=np.zeros((n, 3), dtype=float))
    else:
        colors = data.particles["Color"]

    # Create per-particle radius property if absent
    if "Radius" not in data.particles:
        radii = data.particles_.create_property("Radius", data=np.full(n, 0.18, dtype=float))
    else:
        radii = data.particles["Radius"]

    z = positions[:, 2]

    # Normalize z into [0,1]
    znorm = (z - 0.0) / L
    znorm = np.clip(znorm, 0.0, 1.0)

    # Simple blue -> red gradient by z-position
    # low z = blue, high z = red
    colors[:, 0] = znorm
    colors[:, 1] = 0.2
    colors[:, 2] = 1.0 - znorm

    # Uniform radius for all atoms
    radii[:] = 0.18

pipeline.modifiers.append(PythonScriptModifier(function=setup_particle_visuals))

# -----------------------------
# Set simulation cell explicitly
# -----------------------------
data = pipeline.source.data
data.cell_[:, :] = [
    [L, 0.0, 0.0, 0.0],
    [0.0, L, 0.0, 0.0],
    [0.0, 0.0, L, 0.0],
]

# Turn on PBC flags
data.cell_.pbc = (True, True, True)

# Show simulation cell in rendered output
data.cell.vis.enabled = True
data.cell.vis.rendering_color = (0.1, 0.1, 0.1)

# -----------------------------
# Optional particle display tweaks
# -----------------------------
pipeline.add_to_scene()

scene_data = pipeline.compute()
scene_data.particles.vis.radius = 0.18

# -----------------------------
# Camera/view
# -----------------------------
vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.camera_pos = (1.8 * L, 1.8 * L, 1.2 * L)
vp.camera_dir = (-1.0, -1.0, -0.7)

# -----------------------------
# Render animation
# -----------------------------
vp.render_anim(
    filename=output_file,
    size=(960, 720),
    fps=20,
    background=(1, 1, 1),
    renderer=OpenGLRenderer()
)

print(f"Saved movie to {output_file}")
