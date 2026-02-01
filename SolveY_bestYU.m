function [Y,not_CL] = SolveY_bestYU(P,CL,w,Theorem3)
%% Problem
%  min  yk'Ls yk
%  s.t. yk(jk)=w,yk(ik)=-w,w>0

%%
% Written by YU Wenjun (yuwenjun0203@163.com), written in 2024/3/25, revised in 2024/4/28
%%

if ~exist('w','var')
    w = 1;
end
if ~exist('Theorem3','var')
    Theorem3 = 1;
end

nCL = size(CL,1);
[n,m] = size(P);
Y = zeros(n+m,nCL);

S = sparse(n+m,n+m);
S(1:n,n+1:end) = P;
S(n+1:end,1:n) = P';
%若版本低于2021b版，则使用函数graphconncomp得到最终标签
% [~, y]=graphconyncomp(S);
%若版本高于2021b版，则使用函数conncomp得到最终标签
GG = graph(S);
[y, ~] = conncomp(GG);

not_CL = 0;
for i = 1:nCL
    if y(CL(i,1)) == y(CL(i,2))
        if Theorem3 == 1
            not_CL =  not_CL+1;
            idx_sample = find(y(1:n) == y(CL(i,1)));
            idx_anchor = find(y(n+1:end) == y(CL(i,2)));
            cl = length(find(y(1:CL(i,1)) == y(CL(i,1))));
            cr = length(find(y(1:CL(i,2)) == y(CL(i,2))));
            CL_Mod = [cl,cr];% Modified cannot-link constraint
            P_Mod = P(idx_sample,idx_anchor);
            Y([idx_sample,(n+idx_anchor)],i) = Solve_yk(P_Mod,CL_Mod,w);
        else
            Y(:,i) = Solve_yk(P,CL(i,:),w);
        end
    else
        Y(find(y == y(CL(i,1))),i) = w;
        Y(find(y == y(CL(i,2))),i) = -w;
    end
end
end
