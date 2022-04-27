# Master thesis Marco Codato
In this repository I am going to upload daily the work on my thesis.

## Introduction and overview
I decided to manage data with Python. Since the frames are calibrated and pre-reduced, it is not necessary to use softwares like IRAF.
In the repository are present also the notes that I am writing while doing the work and can be used as a detailed report of the work and an explaination of what you can find in the script source code.

I am going to update my work daily during the week.

## Roadmap

### Done (from the latest)
- I sketched some code that bins each spectrum and compares the frames, chornologically ordered.
- I am starting developing a new script that reads the bkg files generated with the first one.
- **Bkg extraction**. I wrote a script that automatically read and process all the files in a directory.
For each file a preliminary cosmic ray and noise clening is performed. After that all the astronomical sources are detected by looking at the luminosity profile along the slit (wavelength integration).
The remaining areas are saved in a new file, accompained with new information about the procedure in the header.

### To be done soon (in order of priority)
- Develop a strategy to analyze the different spectral features.
- Tune the software using more qualitative criteria

### Next steps (likely):
- Find some references in the literature about the effects of light pollution and in particular concerning spectroscopy.
- ISSUE: if we want to retrieve the sky brightness (insted of the flux only) we have to consider that observations are carried on a portion of the sky delimited by the CCD spatial size (namely the "height")
and the size of the slit. The problem is that the slit size is quite uncertain and heavily affects the precision of the computation.
- Understand how to disentangle the natural sky emission and the contribution from artificial illumination.

### Future major steps (likely):
- Decide how to model the various contributions to the bkg spectrum, in particular LED lights have very different types of spectra.
- Fit the data to obtain the weights of the various sources of the bkg. 
- Develope a convenient strategy to find correlation between the results of the fit and the observation conditions.
- Somehow include error propagation in the whole analysis.