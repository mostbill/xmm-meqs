
# XMM-OM MEQS Project

This code is the technical code content of data reduction, analysis and plot (also potential literature in the future) for the XMM-OM MEQS project.

## Data

- _XMM-OM SUSS_

- _3XMM-DR8_

## Software

- _TOPCAT_

## Usage

- luminosity.py

The entry of the whole code family. Run it to:

1. Group the source from the _fits_ file automatically.

2. Check the _main_ part to:

    - Derive the brevity luminosities, $\alpha_{OX}$ for every source (UV, X-ray).

    - Carry out the linear fitting to the data.

    - Compute the Signal-to-Noise-Ratio (SN ratios) for every source. 

- estimation.py

Sub-code for the _luminosity.py_, as stated in the top of the code.

- wash_fits.py

Sub-code for the _luminosity.py_, as stated in the top of the code.

- match.py

Sub-code for the _luminosity.py_, as stated in the top of the code.