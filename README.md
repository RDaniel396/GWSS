# GWSS
Included in this repo is the files needed to simulate Gravitational Waves to be used as Standard Sirens. The two GWSS_" " are basic python scripts to be added onto the CLASS code pyton wrapper.Both modify the data from CLASS into encorporate the relevant redshift and random fluctuations from perfect LCDM. They also include the error associacted with either detector. Both scripts include a function simulation() this encodes random population of relvant binary mergers that can be measured by the relevant detector. Within the simulation, the relevant signal to noise ratio is determine for each binary. All relevant references are recorded as comments within the script.

# Files
## GWSS_LISA
This script is for LISA. It takes CLASS code, utalises certain columns, modifies the data, includes addtional data to simulate binary mergers. For LISA we have included massive black hole binaries (MBHB) and extream mass ratio inspirals (EMRI). This can be easily updated if relevant ratio/catalogue is provided. This can be chnaged via the mean and standard deviation. 

## GWSS_ET
This script is for the Einstein Telescope (ET). It takes CLASS code, utalises certain columns, modifies the data, includes addtional data to simulate binary mergers. For ET we have included black hole neutron star mergers (BHNS) and binary neutron star (BNS). This can be easily updated if relevant ratio/catalogue is provided.

## lcdm00_background
This file contains the generic background output from CLASS and in this case LambdaCDM. This is just to provide a inital set of data to work with. This will be redundant when .py files are added into the wrapper. 
