clear;clc;
pwdfile = pwd;
addpath(genpath(pwdfile));

datafile = strcat(pwdfile, '\Mvdata&PCs\Mvdata&PCs(ML240_CL80)\');
fileList = dir(fullfile(datafile, '*.mat'));
DataNames = {fileList.name};
ndata=length(DataNames);

fprintf('=============================================\n');
fprintf('            DIRECTORY INFORMATION            \n');
fprintf('=============================================\n');
fprintf('Current working directory:\n');
fprintf('  %s\n\n', pwdfile);
fprintf('Target data directory:\n');
fprintf('  %s\n', datafile);
fprintf('=============================================\n');

fprintf('                    DATASET                  \n');
fprintf('Total files found: %d\n\n', ndata);
for i = 1:ndata
    fprintf('%3d. %s\n', i, DataNames{i});
end
fprintf('\n');
fprintf('=============================================\n');
%%
P = cell(1,ndata);
label = cell(1,ndata);
result = cell(1,ndata);
time = cell(1,ndata);
datanameMatrix = strings(ndata, 1);
nr_train = cell(1,ndata);

data_idx=1:ndata;
Pre = 2;% 样本模一
anchor_rate_list = [0.7,0.6,0.9,0.8];

for i = 1%data_idx
    dataname = DataNames{i}
    datanameMatrix(i,1) = dataname;
%     PCBG(X,c,anchor_rate,ML,CL,gama,Iter,opts)
    load(strcat(datafile,dataname));
    X=M;
    
    %对X进行预处理,样本模一化
    n_View = length(X);
    for v = 1:n_View
        X{v} = Precondition(X{v},Pre);
    end
    
    [n,~] = size(X{1});
    c = length(unique(gnd));
    anchor_rate = anchor_rate_list(i);%锚点率
    gama = 0.1;
    Iter = 30;%外循环的最大迭代数
    opt. style = 0;%DAS锚点选择法
    opt. IterMax = 20;%内循环的最大迭代数
    opt. islocal = 1;%是否局部优化
    opt. Organize_CL = 1;%是否联动must-link和cannot-link约束，整理cannot-link约束
    opt. w = 1;%cnanot-link约束正则项族的超参数
    opt. Theorem2 = 1;%应用Theorem 2进一步优化速度
    opt. Theorem3 = 1;%应用Theorem 3局部更新yk
    opt. Isbreak_Loss = 1;%利用损失函数判断是否达到收敛
    opt. Convergence = 1e-4;
    opt. lambda = 10;
    opt. k = max(3, min(10,fix(anchor_rate*0.5*n/c)));
    %% 检测gnd里面是不是有负数，有就修改
    Min = min(gnd);
    if Min<=0
        gnd = gnd-Min+1;
    end
    %%
    % disp('----------Generation of Paired Constraints----------');
    nML = 0;%要生成的必须连接约束的个数
    nCL = 0;%要生成的不可连接约束的个数
    %[ml, cl] = Generation_of_Paired_Constraints_new(nML, nCL, gnd);

    %% 什么神，什么启动！！！PCBG
    for j = 1:length(ML)
        ml = ML{j};
        cl = CL{j};
        time_start = tic;
        [P{i}, ~,alpha, ~, label{i}, P0, B, Loss,iter] = PCBG(X,c,anchor_rate,ml,cl,gama,Iter,opt);
        
        % ACC, NMI, RandIndx,Purity, Fbeta, Precision, Recall, AdjRandIndx
        [result{i}(j,:)] = ClusteringMeasureTao(gnd, label{i}); 
        result{i}(j,:)
        
        % CPU Time
        time{i}(j,1) = toc(time_start);
        time{i}(j,1)
    
        [nC, nM, rC, rM, rt] = Paired_constraint_validation(label{i}, cl, ml);
        nr_train{i}(j,:) = [nC, nM, rC, rM, rt]; % 实现的CL、ML数量、比例，PCs实现的比例
        nr_train{i}(j,:)
        
        figure;
        imshow(P{i},[]);
        colormap jet;
        colorbar;
    end
    result{i}(j+1,:) = mean(result{i});
    result{i}(j+2,:) = std(result{i}(1:j, :));
    time{i}(j+1,1) = mean(time{i});
    time{i}(j+2,1) = std(time{i}(1:j,1));
    nr_train{i}(j+1,:) = mean(nr_train{i});
    nr_train{i}(j+2,:) = std(nr_train{i}(1:j,:));
end