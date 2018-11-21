# Conventional BCDI (no objective lens involved)

Here are reconstruction scripts for the conventional BCDI experiment without any objective lens. The fringe pattern extends beyond the detector edges, and so we recorded datasets with different (overlapping) detector positions.

Here is a list of the scripts:
- Generate_Mask.m: Construct a mask for the detector, so dead and hot pixels are not used in the reconstruction.
- Reconstruction_Center.m: Reconstruction using a single dataset (the central detector position).
- Processing_Center.m: Correct phase ramp from the diffraction pattern not being centered on the detector, make the real space orthogonal, correct the amplitude for attenuation, calculate the real space and reciprocal space axes.
- Stitch_NoLens.m: Combine the 3 datasets into a single intensity dataset with a higher q-range.
- Reconstruction_Stitched.m: Run the reconstruction on the larger, stitched dataset.
- Processing_Stitched.m: Similar processing to make the result neat, just using a reconstruction from the stitched dataset.
- Simulate_Lenses.m: Simulate the effect of the lens aperture in both the real image and virtual image geometry. The effect is clear as a limit on the q-range in reciprocal space, and a subsequent reduction in the real space resolution.

