function [ML, CL] = Generation_of_Paired_Constraints_new(nML, nCL, gnd)
% Input:
%       - nML: the number of must-link
%       - nCL: the number of cannot-link
%       - gnd: 真实分组

num_cluster = length(unique(gnd));
if nML == 0
    ML = [];
else
    ML = [0, 0];
    nM = 1;
    while nM < nML + 1
        c = randperm(num_cluster, 1);
        logical_index = ismember(gnd, c);
        indices = find(logical_index);
        vector = indices(randperm(length(indices), 2))';
        vector = sort(vector);
        rowIndex = find(ismember(ML, vector, 'rows'));
        if isempty(rowIndex) == 1%说明之前没有生成过这个must-link约束
            nM = nM + 1;
            ML(nM,:) = vector;
        end
    end
    ML = ML(2:nML + 1, :);
end

if nCL == 0
    CL = [];
else
    CL = [0, 0];
    nC = 1;
    while nC < nCL + 1
        c = randperm(num_cluster, 2);
        c1 = c(1);
        c2 = c(2);
        logical_index1 = ismember(gnd, c1);
        indices1 = find(logical_index1);
        randomNumber1 = indices1(randperm(length(indices1), 1));

        logical_index2 = ismember(gnd, c2);
        indices2 = find(logical_index2);
        randomNumber2 = indices2(randperm(length(indices2), 1));

        vector = [randomNumber1, randomNumber2];
        vector = sort(vector);
        rowIndex = find(ismember(CL, vector, 'rows'));
        if isempty(rowIndex) == 1%说明之前没有生成过这个cannot-link约束
            nC = nC + 1;
            CL(nC,:) = vector;
        end
    end
    CL = CL(2:nCL + 1, :);
end
