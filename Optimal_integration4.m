function [y1,y2,P,i,lambda] = Optimal_integration4(B, c, P0, Y, gama, IterMax, islocal)
%% Problem
%  min  ||B - P||^2+gama*Tr(Y'*Ls*Y)+lambda*Tr(F'*LS*F)
%  s.t. P>=0, P*1=1, P'*1=1

% Input:
%       - B: m*n的矩阵，二部图；
%       - c: the number of clusters;
%       - P0: 必须连接约束矩阵；
%       - Y: 保证不可连接约束的指示矩阵，维度为(n+v+m)*p，p表示不可连接约束的个数;
%       - IterMax: the maximum number of iterations;
%       - islocal: 是否局部优化；

% Output:
%       - P: the optimized graph compatible across multiple views;
%       - y1: n个原始样本的标签
%       - y2: m个锚点的标签

zr = 10e-16;
lambda = 10;
[n,m] = size(B);

if ~exist('IterMax','var')
    IterMax=30;
end
if ~exist('islocal','var')
    islocal=1;
end
%%
% disp('----------determine the cluster number----------');%计算合适的聚类数
if ~exist('c','var')
    a = sparse(B);
    a1 = sum(a,2);%按行累加
    D1a = spdiags(1./sqrt(a1),0,n,n);
    a2 = sum(a,1);%按列累加
    D2a = spdiags(1./sqrt(a2'),0,m,m);
    AA1 = D1a*a*D2a;
    
    SS2 = AA1'*AA1;
    SS2 = full(SS2);
    % automatically determine the cluster number
    %%%% 按照特征值的下降趋势来决定聚类数，在下降趋势由快到慢处取聚类数
    [V, ev0, ev]=eig1(SS2,m);   %%% ev包含亲和图的所有特征值，且从大到小排序
    aa = abs(ev); aa(aa>1-zr)=1-eps;
    ad1 = aa(2:end)./aa(1:end-1);   %%%% 后一个特征值比前一个特征值
    ad1(ad1<0.15)=0; ad1 = ad1-eps*(1:m-1)'; ad1(1)=1;   %%% ad1 中值应全部小于等于 1，比值最小处，说明相邻两个特征值差异最大
    ad1 = 1 - ad1; 
    [scores, cs] = sort(ad1,'descend');
    cs = [cs, scores];
    
    c = cs(1);
else
    cs=[];
end
%%
% disp('----------Optimal P----------');%开始计算P
QQ = [P0,B];%初始化Q，n*(m+t)的矩阵
t = size(P0,2);

Dn = sum(QQ,2);%Dn的对角线元素构成的向量，列向量
Dl = spdiags(1./sqrt(Dn),0,n,n);%Dn的(-1/2)次方
Dm = sum(QQ,1);%Dm的对角线元素构成的向量,行向量
Dr = spdiags(1./sqrt(Dm'),0,m+t,m+t);%Dm的(-1/2)次方
X = Dl*QQ*Dr;
X(isnan(X)) = 0;%
X(isinf(X)) = 0;%
XX = X'*X;
XX = full(XX); 
[Fm, ev0, ~]=eig1(XX,c);
Fn=(X*Fm)./(ones(n,1)*sqrt(ev0'));
Fn = (sqrt(2)/2)*Fn;
Fm = (sqrt(2)/2)*Fm;

Yn=Y(1:n,:);% Y的前n行
Ym=Y(end-m+1:end,:);% Y的后m行
disty = L2_distance_1(Yn', Ym');%n*m的距离矩阵

idxa = cell(n,1);
for i=1:n
    if islocal == 1
        idxa0 = find(B(i,:)>0);% 记录B每列里面的非零项的索引
    else
        idxa0 = 1:m;
    end
    idxa{i} = idxa0;
end

idxam = cell(m,1);
for i=1:m
    if islocal == 1
        idxa0 = find(B(:,i)>0);% 记录B每行里面的非零项的索引
    else
        idxa0 = 1:n;
    end
    idxam{i} = idxa0;
end
%%
% disp('----------iterations----------');%开始迭代计算P和F
for i = 1:IterMax
    %%
    %%% update P
    Fn1 = sum(Dl,2).*Fn;
    Fm1 = sum(Dr,2).*Fm;
    distf = L2_distance_1(Fn1',Fm1(end-m+1:end,:)');%n*m的距离矩阵
    Di = gama*disty+lambda*distf;%n*m的距离矩阵
    Q = B-Di/2;
    P = zeros(n,m);%用于按列更新P
    for j = 1:n
        idxa0 = idxa{j}; 
        P(j,idxa0) = EProjSimplex_new(Q(j,idxa0));
    end
    
    P1=zeros(n,m);%用于按行更新P
    for j = 1:m
        idxa0 = idxam{j};
        P1(idxa0,j) = EProjSimplex_new(Q(idxa0,j));
    end
    P=(P+P1)/2;

    %%
    %%% update F
    QQ = [P0,P];
    Dn = sum(QQ,2);%Dn的对角线元素构成的向量，列向量
    newDl = spdiags(1./sqrt(Dn),0,n,n);%Dn的(-1/2)次方
    Dm = sum(QQ,1);%Dm的对角线元素构成的向量,行向量
    newDr = spdiags(1./sqrt(Dm'),0,m+t,m+t);%Dm的(-1/2)次方
    X = newDl*QQ*newDr;
    X(isnan(X)) = 0;%
    X(isinf(X)) = 0;%
    XX = X'*X;
    XX = full(XX);
    
    [newFm, ev0, ev]=eig1(XX,c);
    newFn=(X*newFm)./(ones(n,1)*sqrt(ev0'));
    newFn = (sqrt(2)/2)*newFn;
    newFm = (sqrt(2)/2)*newFm;
    
    %%
    %%% 根据Ls的0特征值数量，改变lambda
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    Dl = newDl;
    Dr = newDr;
    if fn1 < c-zr
        lambda = 2*lambda;
%         Dl = newDl;
%         Dr = newDr;
        Fn = newFn;
        Fm = newFm;
    elseif fn2 > c+1-zr*10
%     elseif fn2 > c+1-0.1
        lambda = 0.8*lambda;
    else
        break;
    end
    
end

%%
% disp('----------Calculate labels----------');%计算标签
%函数graphconncomp在matlab2021b版本中已被移除
% S=sparse(S);
S = sparse(n+m+t,n+m+t);
S(1:n,n+1:end)=[P0,P];
S(n+1:end,1:n)=[P0,P]';

%若版本低于2021b版，则使用函数graphconncomp得到最终标签
% [clusternum, y]=graphconncomp(S);
% y1=y(1:n)';
% y2=y(n+1:end)';

%若版本高于2021b版，则使用函数conncomp得到最终标签
GG = graph(S);
[idx, sizes] = conncomp(GG);
y1=idx(1:n)';
y2=idx(n+1:end)';
end

