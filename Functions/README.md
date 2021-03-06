# Central MatLab functions

This folder contains functions used throughout the other project folders, so I like to keep this folder on the MatLab path so that they are accessible.

Below is a list of the function along with a brief description of their purpose:  
- Be_prop.m: Returns the delta and mu properties of beryllium at the given energy, which must be between 1 and 30 keV.
  - This function loads "Be_delta.txt" and "Be_mu.txt" and performs interpolation.
- bin2.m: Bins a 2D image by the same amount in both directions.
- bin3.m: Bins a 3D image by a specified amount in each dimension, they bin size may vary in different dimensions.
- circ.m: Circle function, assumes a 2D input of the radius, and sets areas of r <= 1 to 1 and 0 otherwise.
- complex2rgb: Converts a 2D complex input into an RGB image, where the amplitude specifies the value and the phase specifies the hue (in HSV).
- CRL_Parameters_1.m: Calculates f, phi, and fN parameters of a CRL.
- CRL_Parameters_2.m: Calculates sigma_D, sigma_a, sigma_v, gamma, and yN parameters for a CRL.
- cycle.m: This is the central 3D BCDI reconstruction function, from the input one can specify which algorithm to use.
  - 1: fienupHIO.m: Simple implementation of HIO, allows a non-binary support (with values in the range [0 1]).
  - 2: error_reduction.m: Simple implementation of ER, allows a non-binary support (with values in the range [0 1]).
  - 3: fienupHIO_2.m: The HIO algorithm, but with contraints on the amplitude for pixels that have been masked. Find a script where cycle(...,3,...) is called to see how it is used.
  - 4: error_reduction_2.m: The ER algorithm, but with contraints on the amplitude for pixels that have been masked. Find a script where cycle(...,4,...) is called to see how it is used.
  - 5: fienupHIO_lm.m: The HIO algorithm with fewer intermediate results, which minimizes the amount of memory needed. However, the error metric can only be calculated in reciprocal space.
  - 6: error_reduction_lm.m: The ER algorithm with fewer intermediate results, which minimizes the amount of memory needed. However, the error metric can only be calculated in reciprocal space.
  - 8: fienupHIOso.m: Original HIO implementation as I received it from Virginie Chamard.
  - 18: error_reduction_so.m: Original ER implementation as I received it from Virginie Chamard.
  - 11: fienupHIO_multi.m: HIO algorithm for a dataset with multiple acquisitions with different masks. The algorithm is run on 1 dataset at a time.
  - 12: error_reduction_multi.m: ER algorithm for a dataset with multiple acquisitions with different masks. The algorithm is run on 1 dataset at a time.
- E2lambda.m: Convert photon energy in eV to wavelength in m.
- edf_info.m: Read header information from .edf files from ID06 ESRF experimental files.
- edf_read.m: Read data from .edf files from ID06 ESRF experimental files.
- edf_write.m: Write data and header in an .edf file format.
- error_metric_data.m: BCDI error metric that compares current guess of the reciprocal space amplitude to the measurement.
- error_metric_real.m: BCDI error metric that compares the real space object between the previous and current iteration.
- error_metric_real.m: BCDI error metric that compares the reciprocal space pattern between the previous and current iteration.
- fconv1.m: 1D convolution implemented with Fourier transforms (fft and ifft).
- fconv2.m: 2D convolution implemented with Fourier transforms (fft2 and ifft2).
- fconvn.m: N-D convolution implemented with Fourier transforms (fftn and ifftn).
- fft3_gpu.m: Implementation of 3D FFT on the GPU for arrays that are potentially larger than the GPU memory. On my configuration this was slower than just running on the CPU, probably due to the data transfer overhead between CPU and GPU memory.
- fileNames.m: Search for file names with a given pattern (e.g. file type) in a specified folder.
- frfft1for.m: 1D (columns) fractional Fourier transform, implemented with a for-loop over the non-transform dimensions. Uses minimal amount of memory, works on the CPU, but is usually slow.
- frfft1gpu.m: 1D (columns) fractional Fourier transform, implemented on the GPU assuming double precision input. Very fast implementation if a GPU is available. Also accepts arrays larger than what can fit on the GPU memory.
- frfft1gpusp.m: 1D (columns) fractional Fourier transform, implemented on the GPU assuming single precision input. Very fast implementation if a GPU is available. Also accepts arrays larger than what can fit on the GPU memory. Single precision computation of FFT is only 2 times faster than double precision, as the computation speed is limited by memory bandwidth.
- frfft1par.m: 1D (columns) fractional Fourier transform, implemented with a parallel for-loop (parfor) over the non-transform dimensions. Uses small amount of memory, works on the CPU, and for small input arrays the overhead is large.
- frfft1vec.m: 1D (columns) fractional Fourier transform, implemented with vectorized code over the non-transform dimensions. Implemented in the CPU, and is the fastest CPU method for small and medium sized arrays. For very large inputs the memory requirements may slow down computation speed.
- frfft2for.m: 2D (columns and rows) fractional Fourier transform. Based on the frfft1for 1D transform.
- frfft2gpu.m: 2D (columns and rows) fractional Fourier transform. Based on the frfft1gpu 1D transform.
- frfft2gpusp.m: 2D (columns and rows) fractional Fourier transform. Based on the frfft1gpusp 1D transform.
- frfft2par.m: 2D (columns and rows) fractional Fourier transform. Based on the frfft1par 1D transform.
- frfft2vec.m: 2D (columns and rows) fractional Fourier transform. Based on the frfft1vec 1D transform.
- FrFT_parameters.m: Given the free space distances and focal lengths of lenses, the fractional Fourier transform parameters may be calculated.
- gauss.m: Calculate a 1D Gaussian with a given FWHM. It is area normalized.
- gaussFWHM.m: Calculate a 1D Gaussian with a given FWHM. Can be area normalzed if specified.
- gauss_RMS.m: Calculate a 1D Gaussian with a given RMS width. Not normalized.
- gpuFFTmaxsize.m: Estimate the computation speed of FFT on the GPU vs array size.
- imshift2.m: Sub-pixel shift of 2D array, using the Fourier shift theorem.
- imshift3.m: Sub-pixel shift of 3D array, using the Fourier shift theorem.
- lambda2E.m: Convert wavelength in m to photon energy in eV.
- load_binary.m: Load a 3D array (real or complex) stored in a binary file. Even though the file size is large due to no compression, reading large files may be faster than the MatLab .mat files. Also see "save_binary.m".
- maskNeighbor2.m: Given a 2D detector input with some masked pixels, for each masked pixel find the mean of the m nearest unmasked pixels.
- maskNeighbor3.m: Given a 3D array input with some masked voxels, for each masked voxel find the mean of the m nearest unmasked voxels.
- optimImagePos.m: In a microscopy setup in imaging geometry, find the exact lens to image plane distance that results in a flat wavefront.
- propAS1.m: 2D wavefront propagation using the angular spectrum method, calls different implementations based on the input array size.
- propAS1cpu.m: 2D wavefront by angular spectrum method implemented on CPU.
- propAS1gpu.m: 2D wavefront by angular spectrum method implemented on GPU. Works with large arrays larger than GPU memory.
- propAS1gpus.m: 2D wavefront by angular spectrum method implemented on GPU. Array must fit on GPU memory.
- propAS1gpuss.m: 2D wavefront by angular spectrum method implemented on GPU. Array must fit on GPU memory.
- propAS.m: 2D wavefront propagation using the angular spectrum method. Original implementation on CPU, allows scaling to vary.
- propFF.m: 2D wavefront propagation using the Fraunhofer approximation. Original CPU implementation.
- propFrFT1.m: 1D wavefront propagation based on the fractional Fourier transform. Core function in the XFrFT package.
- propFrFT2.m: 2D wavefront propagation based on the fractional Fourier transform. Core function in the XFrFT package.
- propIR.m: 2D wavefront propagation based on the impulse response method.
- propTF.m: 2D wavefront propagation based on the transfer function method.
- propTFopt.m: 2D wavefront propagation based on the transfer function method. Optimized implementation on CPU.
- q_range2.m: For a 3D BCDI experiment, calculate the real space and reciprocal space axes. Does not form a full size array, i.e. fast.
- q_range.m: For a 3D BCDI experiment, calculate the real space and reciprocal space axes. Original implementation, forms a full size array, i.e. slow. Prefer "q_range2.m".
- q_vector.m: For a 3D BCDI experiment, calculate the real space and reciprocal space axes. Adapted implementation, but still forms a full size array, i.e. slow. Prefer "q_range2.m".
- rect.m: Rectangle function.
- recta.m: Asymmetric rectangle function.
- save_binary.m: Save a 3D array (real or complex) in a binary file. Even though the file size is large due to no compression, saving large files may be faster than the MatLab .mat files. Also see "load_binary.m".
- shear.m: In 3D BCDI the resulting object is not in an orthogonal frame, and this shear function can transform between measurement space and orthogonal space. The shearing is approximated by integer pixel shifts.
- shear_interp.m: In 3D BCDI the resulting object is not in an orthogonal frame, and this shear function can transform between measurement space and orthogonal space. The shearing is implemented using interpolation.
- shift_amount.m: In 3D BCDI the resulting object is not in an orthogonal frame. This function calculates how much to shear the measurement space to get to an orthogonal frame.
- shrinkwrap.m: Shrinkwrap algorithm for 3D BCDI reconstructions.
- timestr.m: Get the current date and/or time as a string.
- tri.m: Triangle function.
- Vignetting.m: Calculate the CRL RMS vignetting width for use in FrFT propagation, where only an aperture is applied in the object plane (vignetting) and in the lens plane (pupil).
- wavefrontCurvature.m: Calculate the radius of curvature assuming a quadratic phase of the complex input wavefront.
- xcorr2fft.m: 2D cross correlation implemented with Fourier transforms (fft2 and ifft2). Gives the same size output.
- xcorr2fftpad.m: 2D cross correlation implemented with Fourier transforms (fft2 and ifft2). The inputs are padded, and the output can either be padded or same size.



