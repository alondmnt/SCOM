% small scale example script for the SCOM pipeline
% this part consists of building the conservation/divergence networks
% run example.sh for the complete pipeline (including this script)
%
% Alon Diament, Tuller Lab.

load('example/sacCer data');
init_3DND(HiC, bin_index, gene_index, 'example/sacCer 3DND');
init_family_graph('example/sacCer 3DND', 'example/sacCer families');

load('example/schPom data');
init_3DND(HiC, bin_index, gene_index, 'example/schPom 3DND');
init_family_graph('example/schPom 3DND', 'example/schPom families');

build_conserv_net('example/sacCer 3DND fam', 'example/schPom 3DND fam', 15, 'example/conserv net');
build_conserv_net('example/sacCer 3DND fam', 'example/schPom 3DND fam', 25, 'example/conserv net 25');
build_diverg_net('example/sacCer 3DND fam', 'example/schPom 3DND fam', 15, 'example/diverg net');
