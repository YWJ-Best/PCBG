function [Y] = Label_propagation_edit(P,CL,w)
%% Problem
%  min  yk'Ls yk
%  s.t. yk(jk)=w,yk(ik)=-w,w>0

if ~exist('w','var')
    w = 1;
end
nCL = size(CL,1);
[n,m]=size(P);
Y = zeros(n+m,nCL);
yu = [1;-1];
dn = sum(P,2);
dm = sum(P',2);
inv_dn = dn.^(-1);
Dm = diag(dm);

for i = 1:nCL
    C = CL(i,:);
    Max=max(C);
    Min=min(C);
    idx_i = setdiff(1:n, C);
    P1 = P([Min,Max],:);
    P2 = P(idx_i,:);

    inv_dn_2 = inv_dn(idx_i);
    inv_Dn_2_P2 = P2.*inv_dn_2;
    Q = Dm-P2'*inv_Dn_2_P2;
    beta = 1e-6;
    inv_Q = inv(Q+beta*eye(m));
    B = inv_Q*P1';
    A = inv_Dn_2_P2*B;

    yd = [A;B]*yu;
    Y(:,i)=[yd(1:Min-1);w;yd(Min:Max-2);-w;yd(Max-1:end)];

    % idx_i_all = setdiff(1:(n+m), C);
    % Y(Min,i) = 1;
    % Y(Max,i) = -1;
    % Y(idx_i_all,i) = [A;B]*yu;
    
%     Y(Y(:,i) > 1,i) = 1;
%     Y(Y(:,i) < -1,i) = -1;
    
end
end
