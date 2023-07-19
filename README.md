# User Manual: Downward Continuation and Redatuming to the Seafloor of Streamer Data

We present a software for the Downward Continuation (DC) of marine field data. The DC is a virtual transformation of the shot gathers from datum 1 (shots and receivers at the sea surface) to datum 2 (shots and receiver at seafloor). This redatuming enhances the refractions present in the data, which are no longer hidden or overlapped by the seafloor reflections. The DC is done in two steps:

- In the first step or DC<sub>1</sub>, the receivers at each shot gather are downward continued to datum 2.
- In the second step or DC<sub>2</sub>, the shot gathers are downward continued to datum 2.

## Inputs

As inputs, the following are required:
- Shot gathers in SU format.
- Bathymetry at each shot gather in an ASCII file.
- A series of parameters provided through an input file, which is read during execution.
- By default, the p-wave velocity model for the water column is considered constant. However, XBT data can be provided to build a specific Vp model.

For more details on the physics implemented in the software, see the reference pre-print [1]. This reference is also revealing to know which type of streamer data is a good candidate for downward continuation, and also to understand the validity limits of the DC results.

## Installation

Instructions on how to install and use this software are available in the PDF manual located in this directory and named: [UserGuide_DC.pdf](UserGuide_DC.pdf)

## Development

- Development is hosted on the GitHub repository: [github/ejimeneztejero/DC](https://github.com/ejimeneztejero/DC).

## Author

- The author is Clara Estela Jiménez Tejero.
- This software was developed at Barcelona Center for Subsurface Imaging, at ICM-CSIC.

## Reference

Clara Estela Jimenez-Tejero et al, Downward continuation of marine seismic reflection data: an undervalued tool to improve velocity models, Geophysical Journal International, Volume 230, Issue 2, August 2022, Pages 831–848, https://doi.org/10.1093/gji/ggac087
