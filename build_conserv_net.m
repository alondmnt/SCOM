function build_conserv_net(infile1, infile2, k, outfile)
% build_conserv_net(infile1, infile2, k, outfile)
%   generate a conservation network based on common top/bottom-k (percent) 3DND.
% INTPUT:
%   [infile1] path to family graph of species-1.
%   [infile2] path to family graph of species-2.
%   [k] percentile threshold.
%   [outfile] path to file.
% OUTPUT:
%   [G] conservation network (green edges)
%   [R] conservation network (red edges)
%
% Alon Diament, Tuller Lab.

D1 = load(infile1, 'D');
D2 = load(infile2, 'D');

coords = find(any(D1.D, 1) & any(D2.D, 1));  % common families
D1 = full(D1.D(coords, coords));
D2 = full(D2.D(coords, coords));
D1(D1 == 0) = NaN;  % missing values
D2(D2 == 0) = NaN;

[testnan, isort1] = sort(D1(:), 'ascend');
isort1 = isort1(~isnan(testnan));
[testnan, isort2] = sort(D2(:), 'ascend');
isort2 = isort2(~isnan(testnan));

k = ceil(k/100 * length(isort1));

fprintf('conserved CLOSE (CoLoc/green edges)\n');
% select values that are close in both species
ind = intersect(isort1(1:k), isort2(1:k));
G = false(size(D1));
G(ind) = 1;
G = (G + G') > 0;  % enforce symmetry

fprintf('conserved FAR (CoFar/red edges)\n');
% select values that are far in both species
ind = intersect(isort1(end-k+1:end), isort2(end-k+1:end));
R = false(size(D1));
R(ind) = 1;
R = (R + R') > 0;

save(outfile, 'G', 'R', 'coords', '-v7');
fprintf('saved %s\n', outfile);
