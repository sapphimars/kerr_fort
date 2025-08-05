import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Load data
with open("image.dat", "r") as f:
    lines = f.readlines()

data = []
for line in lines[1:]:
    if line.strip():
        row = [float(x) for x in line.split()]  # Note: float now
        data.extend(row)

image = np.array(data).reshape((256, 256))

# Create enhanced visualization
fig, ax1 = plt.subplots(figsize=(10, 10))

# Main image with enhanced contrast
im1 = ax1.imshow(
    image,
    cmap="hot",
    origin="lower",
    vmin=0,
    # vmax=np.percentile(image[image > 0], 95),
)
threshold = 0.01
overlay_mask = (image < threshold) & (image > 0)
red_overlay_data = np.full(image.shape, np.nan)
red_overlay_data[overlay_mask] = 1
cmap = colors.ListedColormap(["white"])

# im2 = ax1.imshow(red_overlay_data, origin="lower", vmin=0, cmap=cmap, alpha=1)


ax1.set_title("Kerr Black Hole (a=0.0)")
ax1.set_xlabel("Pixel X")
ax1.set_ylabel("Pixel Y")
plt.colorbar(im1, ax=ax1, label="Intensity")
plt.show()
"""
fig, ax2 = plt.subplots(figsize=(10, 10))
# Contour plot to show structure
im2 = ax2.contourf(
    image,
    levels=100,
    cmap="plasma",
    vmin=0,
    vmax=np.percentile(image[image > 0], 95),
    origin="lower",
)
ax2.set_title("Contour View - Structure Analysis")
ax2.set_xlabel("Pixel X")
ax2.set_ylabel("Pixel Y")
plt.colorbar(im2, ax=ax2, label="Intensity")

# Add annotations for expected features
ax1.text(
    0.02,
    0.98,
    "Expected:\n• Central shadow\n• Photon sphere\n• Primary disk\n• Secondary image",
    transform=ax1.transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
)

# plt.tight_layout()
plt.show()
"""
