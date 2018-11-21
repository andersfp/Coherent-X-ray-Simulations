# Objective BCDI (virtual image geometry)

Here are many files, as I have made many different attempts to combine different datasets (most were not successful). We acquired 9 datasets, moving the objective lens (CRL, N = 20) in a square 3x3 grid. Below I will only describe the most important scripts.

Reconstruction scripts:
- Reconstruction_Center.m: A normal reconstruction using a single dataset.
- Reconstruction_Ptycho_1.m: Fourier synthesis reconstruction named "parallel" in [1].
- Reconstruction_Ptycho_5.m: Fourier synthesis reconstruction named "serial" in [1].

The other Reconstruction_Ptycho_N.m scripts are different ideas of obtaining a combined reconstruction, but they are not successful. For Reconstruction_Stiched_N.m I tried to pre-stitch the data into a single dataset (the stitching is in Stitch_Virtual_N.m), but that reconstruction does not work either. I also tried aligning the data by shifting the datasets to the correct position before the reconstruction (Shift_Virtual.m and Reconstruction_Shifted.m), but this was not successful either (one potential problem was that the shifting was only done in integer number of pixels, which is not accurate enough).

In order to get the real space object with correct axis, I used the following scripts:
- Find_Magnification_2.m: Compare the detector fringe pattern to that of the conventional BCDI experiment to determine the magnification of the fringe pattern.
- Lens_Parameters.nb: A Mathematica notebook for calculating the lens position from the detector fringe magnification.
- Processing_X_N.m: Correct phase ramp from the diffraction pattern not being centered on the detector, make the real space orthogonal, correct the amplitude for attenuation, calculate the real space and reciprocal space axes.

I did addtional testing of the MTF, shifting data using the Fourier shifting theorem, and a bunch of other experiments.

The two PowerPoints contain results on a few of these tests.

The fitting results Virtual_Center_ERF_Fit.sfit can be opened with the 'cftool' (built-in MatLab curve fitting UI tool). Here I have fitted the top edge of the object by an error function to estimate the resolution. It contains data from reconstruction without objective lens and for the real image geometry as well.

[1]  
"X-ray coherent diffraction imaging with an objective lens: towards 3D mapping of thick polycrystals".  
A. F. Pedersen, V. Chamard, C. Detlefs, T. Zhou, D. Carbone, H. F. Poulsen.  
arXiv:1810.04268 [physics.ins-det]
