import os
from PIL import Image

# Path to the folder with PNGs
folder_path = 'frames_for_gif_real_tles_only'

# Get all .png files, sorted by filename
frames = sorted(
    [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.png')]
)

# Load with Pillow
frames_pil = [Image.open(f) for f in frames]

# Save as GIF with custom durations
frames_pil[0].save(
    'heatmaps_real_tles.gif',
    save_all=True,
    append_images=frames_pil[1:] + [frames_pil[-1]],  # repeat last frame
    duration=[250] * len(frames_pil) + [2500],        # ms per frame
    loop=0
)
