# BCDI through the CRL objective pass-through

These reconstructions are just sanity checks that the reconstructions still works when the objective lens box is present (which gives more air and some kapton foil in the beam).

The mask files are binary image files showing which portions of the detector are usable (not blocked by the objective lens). The MaxiPix mask is used to filter out dead and hot pixels on the detector.

Reconstruction_N runs the BCDI reconstruction with a CRL lens box with either 20 or 70 lenslets. These should give the same result, as the lens box is the same.

