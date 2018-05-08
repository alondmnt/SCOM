#!/bin/bash

### NOTE: for the sake of the example, we set min_seed to 0 in scom.py ###
### the original value in the paper (for a complete dataset) was 20    ###
### also, here we do not generate an empirical distribution of edges   ###

# prepare conservation / divergence networks
matlab -nodesktop -nodisplay -r "example;quit"

# conserved SCOM
python scom.py con "example/conserv net"
# divergent SCOM in org1
python scom.py div "example/diverg net 1"
# divergent SCOM in org2
python scom.py div "example/diverg net 2"
# separated SCOM
python scom.py sep "example/conserv net 25"

