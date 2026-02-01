function [ML, CL] = Generation_of_Paired_Constraints(nML, nCL, gnd)
% Input:
%       - nML: the number of must-link
%       - nCL: the number of cannot-link
%       - gnd: 真实分组

n=size(gnd,1);
if nML == 0
    ML = [];
    if nCL == 0
        CL = [];
    else
        CL = [0,0];
        nC=1;
        while nC<nCL+1
            randomNumber1 = randi([1, n-1]);
            randomNumber2 = randi([randomNumber1+1, n]);
            vector=[randomNumber1, randomNumber2];
            if gnd(randomNumber1) ~= gnd(randomNumber2)
                rowIndex = find(ismember(CL, vector, 'rows'));
                if isempty(rowIndex) == 1%说明之前没有生成过这个cannot-link约束
                    nC=nC+1;
                    CL(nC,:)=vector;
                end
            end
        end
        CL = CL(2:nCL+1,:); 
    end
else
    if nCL == 0
        CL = []; 
        ML = [0,0];
        nM=1;
        while nM<nML+1
            randomNumber1 = randi([1, n-1]);
            randomNumber2 = randi([randomNumber1+1, n]);
            vector=[randomNumber1, randomNumber2];
            if gnd(randomNumber1) == gnd(randomNumber2)
                rowIndex = find(ismember(ML, vector, 'rows'));
                if isempty(rowIndex)==1%说明之前没有生成过这个must-link约束
                    nM=nM+1;
                    ML(nM,:)=vector;
                end
            end
        end
        ML=ML(2:nML+1,:);
    else
        ML=[0,0];
        CL=[0,0];
        nM=1;
        nC=1;
        while nM<nML+1 || nC<nCL+1
            randomNumber1 = randi([1, n-1]);
            randomNumber2 = randi([randomNumber1+1, n]);
            vector=[randomNumber1, randomNumber2];
            if gnd(randomNumber1) == gnd(randomNumber2)
                rowIndex = find(ismember(ML, vector, 'rows'));
                if isempty(rowIndex)==1%说明之前没有生成过这个must-link约束
                    nM=nM+1;
                    ML(nM,:)=vector;
                end
            else
                rowIndex = find(ismember(CL, vector, 'rows'));
                if isempty(rowIndex)==1%说明之前没有生成过这个cannot-link约束
                    nC=nC+1;
                    CL(nC,:)=vector;
                end
            end   
        end
        ML=ML(2:nML+1,:);
        CL=CL(2:nCL+1,:);
    end
end

end

