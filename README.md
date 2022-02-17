# User Manual: Downward Continuation and redatuming to the seafloor of streamer data

This is the first version of the software for the Downward Continuation (DC) of marine field data. The DC is a virtual transformation of the shot gathers from datum 1 (shots and receivers at the sea surface) to datum 2 (shots and receiver at seafloor).
This redatuming enhance the refractions present in the data, which are not anymore hidden or overlapped by the seafloor reflections. The DC is done in two steps:

- In the first step or DC<sub>1</sub>, the receivers at each shot gather are downward continued to datum 2.
- In the second step or DC<sub>2</sub>, the shot gathers are downward continued to datum 2.

As inputs:
- The shot gathers are required in SU format.
- The bathymetry at each shot gather in a ascii file.
- Also, a series of parameters have to be provided through an input file which is read in the execution line.
- By default, the p-wave velocity model for the water column is considered constant but also XBT data can be provided to build a specific Vp model.

For more details on the physics implemented in the software, see the reference pre-print [1]. This reference is also revealing to know which type of streamer data is a good candidate for downward continuation, and also to understand the validity limits of the DC results.

## Installation
Instructions on how to install and use this software are available in the PDF manual located in this directory and named: [UserGuide.pdf](UserGuide.pdf)

## Development
- Development is hosted on GitHub repository:
[github/ejimeneztejero/DC](https://github.com/ejimeneztejero/DC).

## Author
- The author is Clara Estela Jiménez Tejero.
- This software was developed at Barcelona Center for Subsurface Imaging, at ICM-CSIC.

## Reference
- Clara Estela Jimenez Tejero, Cesar R. Ranero, Valenti Sallares and Claudia Gras. 'Open source downward continuation to the seafloor of streamer data', arXiv: https://arxiv.org/abs/2106.00646.
