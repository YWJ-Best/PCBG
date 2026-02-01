function [ML, ML_group, t, CL, CL_group, CL_group1] = Organize_PCsYU(ML, CL)
% Inpu
%       - ML: 规划/整理 之前的must-link
%       - CL: cannot-link

% Output:
%       - ML: 规划/整理 之后的must-link，此时的must-link约束与原本的must-link约束完全等价，但数目会更少
%       - ML_g: 储存了must-link指出的每个族群中所有样本的序号
%       - t: must-link指出的族群数量
%       - CL: 联合must-link约束规划/整理 之后的cannot-link，
%             此时的cannot-link约束与原本的cannot-link约束在Subgroup-level等价，但数目会更少
%       - CL_group: cannot-link指出的MLc中的敌对子群的序号对，即MLc中子群的不可连接约束
%       - CL_group1: 对于cannot-link中只有一个对象落在sub-group中时，将其转化为样本-sub-group的cannot-lnik，
%                    第一列为子群序号，第二列为样本序号

%   Written by Yu Wenjun (YUWENJUN0203@163.com), written in 2023/12/12, revised in 2024/5/20

nML = size(ML,1);
nCL = size(CL,1);
ML_group = cell(1,1);
a = 0;%开启总循环，说明ML还有多余的，没有被完全 规划/整理
M=[];
t = 0;

if nML > 0 % 有must-link约束
    while a == 0
        ML_r = ML(1,:);
        ML(1,:) = [];
        v = ML_r;%用于储存，must-link指出的一个连通子图内的所有索引
        b = 0;
        while b ==0
            b = 1;
            ML_r1 = [];%用于储存，与ML_r中样本位于同一个连通子图的样本的索引
            idxa = [];
            for j = 1:size(ML,1)
                if length(setdiff(ML(j,:), ML_r)) ~=2
                    ML_r1 = [ML_r1, ML(j,:)];
                    idxa = [idxa, j];
                    b = 0;
                end
            end
            v = [v, ML_r1];
            ML_r = unique(ML_r1);
            if size(idxa,1) ~= 0 
                ML(idxa,:) = [];
            end
        end
        if size(ML,1) == 0
            a = 1;
        end
        v = unique(v);%去重并从小到大排序
        t = t+1;
        ML_group{t} = v;
        % 现通过该连通子图内的所有索引v，重新整理形成新的must-link
        k = 0;
        for m = v(1:end-1)
            k = k+1;
            M = [M; [v(k),v(k+1)]];
        end
    end
end
ML = M;

insubgroup_idx = [];%cannot-link约束中两个样本都落在某个subgroup
insubgroup_idx1 = [];%cannot-link约束中只有一个样本落在某个subgroup
if nML > 0 % 有must-link约束
    ML_samples = unique([ML(:,1),ML(:,2)]);
    for i = 1:nCL
        Max = max(CL(i,:));
        Min = min(CL(i,:));
        isme = sum(ismember(ML_samples, Min)) + sum(ismember(ML_samples, Max));
        if isme == 2
            insubgroup_idx = [insubgroup_idx,i];
        elseif isme == 1
            insubgroup_idx1 = [insubgroup_idx1,i];
        end
    end
end
n_insub = size(insubgroup_idx,2);
move_idx = [];
CL_group = [];
if n_insub >= 2
    Max = max(CL(insubgroup_idx(1),:));
    Min = min(CL(insubgroup_idx(1),:));
    for j = 1:t
        if sum(ismember(ML_group{j}, Min)) == 1
            CL_group(1,1) = j;
        elseif sum(ismember(ML_group{j}, Max)) == 1
            CL_group(1,2) = j;
        end
    end
    for i = 2:n_insub
        Max = max(CL(insubgroup_idx(i),:));
        Min = min(CL(insubgroup_idx(i),:));
        for j = 1:t
            if sum(ismember(ML_group{j}, Min)) == 1
                CL_group(i,1) = j;
            elseif sum(ismember(ML_group{j}, Max)) == 1
                CL_group(i,2) = j;
            end
        end
        if sum(ismember(CL_group(1:i-1,:), CL_group(i,:), 'rows')) >= 1 %说明前面出现了子群级的等价cannot-link
            move_idx = [move_idx,insubgroup_idx(i)];
        end
    end
end

n_insub1 = size(insubgroup_idx1,2);
CL_group1 = [];
if n_insub1 >= 2% if n_insub1 <=1, we don't have to do anything
    Max = max(CL(insubgroup_idx1(1),:));
    Min = min(CL(insubgroup_idx1(1),:));
    for j = 1:t
        if sum(ismember(ML_group{j}, Min)) == 1
            CL_group1(1,1) = j;
            CL_group1(1,2) = Max;
        elseif sum(ismember(ML_group{j}, Max)) == 1
            CL_group1(1,1) = j;
            CL_group1(1,2) = Min;
        end
    end
    for i = 2:n_insub1
        Max = max(CL(insubgroup_idx1(i),:));
        Min = min(CL(insubgroup_idx1(i),:));
        for j = 1:t
            if sum(ismember(ML_group{j}, Min)) == 1
                CL_group1(i,1) = j;
                CL_group1(i,2) = Max;
            elseif sum(ismember(ML_group{j}, Max)) == 1
                CL_group1(i,1) = j;
                CL_group1(i,2) = Min;
            end
        end
        if sum(ismember(CL_group1(1:i-1,:), CL_group1(i,:), 'rows')) >= 1 %说明前面出现了子群级的等价cannot-link
            move_idx = [move_idx,insubgroup_idx1(i)];
        end
    end
end
CL(move_idx,:) = [];
end

