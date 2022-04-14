# Master thesis Marco Codato

## Introduction and overview
In this repository I am going to upload daily the work on my thesis.

I decided to manage data with Python. Since the frames are calibrated and pre-reduced, it is not necessary to use softwares like IRAF.
In the repository are present also the notes that I am writing while doing the work and can be used as a detailed report of the work and an explaination of what you can find in the script source code.

## Roadmap

### Done (from the latest)
- A quantitative approach in the source selection is quite complicated due to the possible luminosity profiles of different sources along the slit. It is more practical to use the current empirical approach.
- In the header of the new file is printed the time when the file was processed. Other information like the extracted rows is unecessary since the shape of the new file is the same of the original one and the masked regions are simply mixels with no numerical values.
- The script save the bkg masked spectrum in a mew file, mantaining the orginal header too.
- Write a script that opens a FITS file, automatically detects and remove astornomical sources.

### To be done soon (in order of priority)
- I need more images to test the script and check its robustness.

### Next steps (likely):
- Find some references in the literature about the effects of light pollution and in particular concerning spectroscopy.
- Test the script with more frames to check the robustness.
- Decide wether integrate the bkg over the slit position or not.
- ISSUE: if we want to retrieve the sky brightness (insted of the flux only) we have to consider that observations are carried on a portion of the sky delimited by the CCD spatial size (namely the \height")
and the size of the slit. The problem is that the slit size is quite uncertain and heavily affects the precision of the computation.
- Understand how to disentangle the natural sky emission and the contribution from artificial illumination.

### Future major steps (likely):
- Decide how to model the various contributions to the bkg spectrum, in particular LED lights have very different types of spectra.
- Fit the data to obtain the weights of the various sources of the bkg. 
- Develope a convenient strategy to find correlation between the results of the fit and the observation conditions.