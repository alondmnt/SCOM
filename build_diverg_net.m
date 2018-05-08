function build_diverg_net(infile1, infile2, k, outfile)
% build_diverg_net(infile1, infile2, k, outfile)
%   generate a divergence network based on common bottom-k (percent) 3DND.
% INTPUT:
%   [infile1] path to family graph of species-1.
%   [infile2] path to family graph of species-2.
%   [k] percentile threshold.
%   [outfile] path to file.
% OUTPUT:
%   [G] side-1 divergence network, written to [outfile] with suffix 1.
%   [G] side-2 divergence network, written to [outfile] with suffix 2.
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

k = ceil(k/100 * numel(D1));

fprintf('\ndiverged edges, co-loc in organism (1)\n');
ind = intersect(isort1(1:k), isort2(end-k+1:end));
G = false(size(D1));
G(ind) = 1;
G = (G + G') > 0;
out1 = [outfile, ' 1'];
save(out1, 'G', 'coords', '-v7');
fprintf('saved %s\n', out1);

fprintf('diverged edges, co-loc in organism (2)\n');
ind = intersect(isort2(1:k), isort1(end-k+1:end));
G = false(size(D1));
G(ind) = 1;
G = (G + G') > 0;
out2 = [outfile, ' 2'];
save(out2, 'G', 'coords', '-v7');
fprintf('saved %s\n', out2);
