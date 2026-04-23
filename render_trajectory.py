from ovito.io import import_file
from ovito.vis import Viewport, OpenGLRenderer

# ---- User settings ----
input_file = "trajectory.xyz"
output_file = "argon_movie.mp4"

# Box length from your Rahman-style run
L = 10.229343992174416

# ---- Load trajectory ----
pipeline = import_file(input_file, multiple_frames=True)

# Set simulation cell explicitly
data = pipeline.source.data
data.cell_[:, :] = [
    [L, 0.0, 0.0, 0.0],
    [0.0, L, 0.0, 0.0],
    [0.0, 0.0, L, 0.0],
]

# ---- Set up camera/view ----
vp = Viewport()
vp.type = Viewport.Type.Perspective

# Put camera outside the box looking toward center
vp.camera_pos = (1.8 * L, 1.8 * L, 1.2 * L)
vp.camera_dir = (-1.0, -1.0, -0.7)

# ---- Render movie ----
vp.render_anim(
    filename=output_file,
    size=(960, 720),
    fps=20,
    background=(1, 1, 1),
    renderer=OpenGLRenderer()
)

print(f"Saved movie to {output_file}")
