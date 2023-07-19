# User Manual: Downward Continuation and Redatuming to the Seafloor of Streamer Data

This README provides an overview and instructions for using the Downward Continuation (DC) software for marine field data. DC is a virtual transformation of shot gathers from datum 1 (shots and receivers at the sea surface) to datum 2 (shots and receivers at the seafloor). This redatuming process enhances the visibility of refractions in the data by removing the interference caused by seafloor reflections. The DC process consists of two steps:

- DC<sub>1</sub>: Downward continuation of the receivers at each shot gather to datum 2.
- DC<sub>2</sub>: Downward continuation of the shot gathers to datum 2.

## Inputs

To use the software, you need the following inputs:

- Shot gathers in SU format.
- Bathymetry data for each shot gather in an ASCII file.
- An input file containing various parameters that are read during execution.
- By default, the software assumes a constant p-wave velocity model for the water column. However, you can also provide XBT data to create a specific Vp model.
- For a more detailed understanding of the underlying physics and to determine the suitability and limitations of DC results for different streamer data, refer to the reference [1]


## Installation

- To install and use the software, please refer to the instructions provided in the PDF manual located in this directory:  [UserGuide_DC.pdf](UserGuide_DC.pdf)

## Development

- The development of this software is hosted on the GitHub repository: 
  [github/ejimeneztejero/DC](https://github.com/ejimeneztejero/DC).

## Author

- Clara Estela Jiménez Tejero is the author of this software.
- The software was developed at the Barcelona Center for Subsurface Imaging, ICM-CSIC.

## Reference

[1] Clara Estela Jimenez-Tejero et al, Downward continuation of marine seismic reflection data: an undervalued tool to improve velocity models, Geophysical Journal International, Volume 230, Issue 2, August 2022, Pages 831–848, https://doi.org/10.1093/gji/ggac087
