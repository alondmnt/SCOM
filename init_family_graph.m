function init_family_graph(infile, family_set)
% init_family_graph(infile, family_set)
%   here we take a variable (distances) in gene coordinates and transform 
%   it to (universal) families coordinates, where each family is a node.
%
% INPUT:
%   [infile]: a file containing the variable [D].
%   [family_set]: a file containing the families dataset, i.e., a Nx1 cell array
%       [orto] with a list of gene indices in each cell/family (can be empty). 
%
% OUPUT: (saved to file with a 'fam' suffix)
%   [D]: a NxN matrix of distances between families.
%
% Alon Diament, Tuller Lab.

load(family_set, 'orto');
fprintf('\nloaded family_set %s\n', family_set);
load(infile);
fprintf('loaded distances\n');

DataFile.VarName = 'D';
DataFile.FileName = infile;
DataFile.FamilySet = family_set;
DataFile.Function = 'init_family_graph';

[outfile{1}, outfile{2}, outfile{3}] = fileparts(DataFile.FileName);

outfile = fullfile(outfile{1}, [outfile{2}, ' fam', outfile{3}]);

n_node = length(orto);
n_inter = sum(~cellfun(@isempty, orto));
n_inter = ceil(n_inter^2 / 2);
inters = zeros(n_inter, 3);
idx = 1;
% generating a sparse table of average interactions between families
for i = 1:n_node
    i_set = orto{i};
    if isempty(i_set)
        continue
    end
    Di = D(i_set, :);
    inner_calc = zeros(i, 3);
    for j = 1:i  % can use parfor
        j_set = orto{j};
        if isempty(j_set)
            continue
        end
        d = Di(:, j_set);
        div = nnz(~isnan(d) & d~=0);  % unknown values are NaN/zero and ignored
        if div == 0
            continue
        end
        inner_calc(j, :) = [i, j, nansum(d(:)) / div];
    end
    
    inner_calc = inner_calc(inner_calc(:,3) ~= 0, :);
    len = size(inner_calc, 1);
    inters(idx:idx + len - 1, :) = inner_calc;
    idx = idx + len;
end

inters = inters(inters(:, 3) ~= 0, :);
D = sparse(inters(:, 1), inters(:, 2), inters(:, 3), ...
           n_node, n_node);
D = D + D' - diag(diag(D));
fprintf('buily family graph\n');

save(outfile, 'D', 'DataFile', '-v7.3');
fprintf('saved %s\n', outfile);
