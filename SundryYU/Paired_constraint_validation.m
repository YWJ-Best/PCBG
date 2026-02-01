function [nC, nM, rC, rM, rt, CL_un, ML_un, CL_yes, ML_yes] = Paired_constraint_validation(label, CL, ML)


% output:
%         Cl:没被完成的CL约束
%         Ml:没被完成的ML约束

nCL = size(CL,1);
nML = size(ML,1);

nC = 0;
nM = 0;
CL_un=[];
ML_un=[];
CL_yes=[];
ML_yes=[];
if nCL>0
    for i = 1:nCL
        if label(CL(i,1)) ~= label(CL(i,2))
            nC = nC+1;
            CL_yes = [CL_yes;CL(i,:)];
        else
            CL_un = [CL_un;CL(i,:)];
        end
    end
    rC = nC/nCL;
else
    rC="无cannot-link";
end
if nML>0
    for i = 1:nML
        if label(ML(i,1)) == label(ML(i,2))
            nM = nM+1;
            ML_yes = [ML_yes;ML(i,:)];
        else
            ML_un = [ML_un;ML(i,:)];
        end
    end
    rM = nM/nML;
else
    rM="无must-link";
end
if nCL+nML>0
    rt = (nC+nM)/(nCL+nML);
else
    rt = "无成对约束";
end
