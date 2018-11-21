# Objective BCDI (real image geometry)

This is a reconstruction with the real image geomtry, which significantly limits the resolution. The objective lens was a CRL with N = 70 lenslets. We recorded several datasets with different objective lens postions. However, we only got around to do the reconstruction on a single dataset (see the virtual geometry folder for reconstructions with multiple datasets). The basic functionality is:
- Reconstruction_Center.m: Perform a reconstruction using a single dataset.
- Processing_Center.m: Correct phase ramp from the diffraction pattern not being centered on the detector, make the real space orthogonal, correct the amplitude for attenuation, calculate the real space and reciprocal space axes.

Additional parameters needed for getting a correct real space axis:
- Find_Magnification_2.m: Compare the detector fringe pattern to that of the conventional BCDI experiment to determine the magnification of the fringe pattern (version 2 of the script should give superior accuracy).
- Lens_Parameters.nb: A Mathematica notebook for calculating the lens position from the detector fringe magnification.

