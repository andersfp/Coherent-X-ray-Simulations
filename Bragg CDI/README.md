# Simulations for Objective BCDI

Before the objective BCDI experiment we did some simulations to verify how it works. We simulate the forward propagation with and without a CRL objective lens (Lens vs Free_Space), and in each case we consider two different noise levels (1 vs 2). The simulations are documentet in the journal article:  
"Numerical study of Bragg CDI on thick polycrystalline specimens".  
A. F. Pedersen, V. Chamard, and H. F. Poulsen.  
Optics Express, vol. 26, No. 18 (2018), p. 23411-23425.  
DOI: https://doi.org/10.1364/OE.26.023411

The main simulations can all be run through the Run_All.m script, and they are as follows:
- Experimental_Parameters.m: Definition of all the distances, lens parameters, number of pixels, etc.
- Generate_Object.m: Make the phantom object for the simulations.
- Generate_Exit_Field.m: Calculate the 2D complex wavefield at the object plane for each tilt angle of the sample.
- Forward_Projection.m: Propagate the exit fields to the detector, both with and without an objective lens. The detector intensity gets noise added, using two different signal levels.
- Lowres_Shrinkwrap_Free_Space_1.m: Perform a low resolution reconstruction of the simulation without an objective lens with high signal level. The low resolution reconstruction is performed on the GPU.
- Lowres_Shrinkwrap_Free_Space_2.m: Perform a low resolution reconstruction of the simulation without an objective lens with low signal level. The low resolution reconstruction is performed on the GPU.
- Reconstruction_Shrinkwrap_Free_Space_1.m: Perform the full resolution reconstruction of the simulation without an objective lens with high signal level. The reconstruction is performed on the CPU, but using the low resolution result from the GPU as a seed.
- Reconstruction_Shrinkwrap_Free_Space_2.m: Perform the full resolution reconstruction of the simulation without an objective lens with low signal level. The reconstruction is performed on the CPU, but using the low resolution result from the GPU as a seed.
- Lowres_Shrinkwrap_Lens_1.m: Perform a low resolution reconstruction of the simulation with a CRL objective lens with high signal level. The low resolution reconstruction is performed on the GPU.
- Lowres_Shrinkwrap_Lens_2.m: Perform a low resolution reconstruction of the simulation with a CRL objective lens with low signal level. The low resolution reconstruction is performed on the GPU.
- Reconstruction_Shrinkwrap_Lens_1.m: Perform the full resolution reconstruction of the simulation with a CRL objective lens with high signal level. The reconstruction is performed on the CPU, but using the low resolution result from the GPU as a seed.
- Reconstruction_Shrinkwrap_Lens_2.m: Perform the full resolution reconstruction of the simulation with a CRL objective lens with low signal level. The reconstruction is performed on the CPU, but using the low resolution result from the GPU as a seed.

The scripts N_Fig_x.m produce the figures found in the article above.

The remaining scripts are used to estimate resolution and to perform various sanity checks.

