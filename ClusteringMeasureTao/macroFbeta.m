function [ macrofbeta ] = macroFbeta( y_pre, y_true, beta )
% computer the macro fbeta score
%  
if ~exist('beta','var')
    beta = 1;
end
y_true = y_true(:);
y_pre = y_pre(:);
mat=confusionmat(y_true, y_pre);
classNum = size(mat,1);
Fbeta = zeros(classNum,1);
P = zeros(classNum,1);
R = zeros(classNum,1);
for i = 1:classNum
    P(i) = mat(i,i)/sum(mat(:,i));
    R(i) = mat(i,i)/sum(mat(i,:));
    if P(i) == 0 || R(i) == 0
        Fbeta(i) = 0;
    else
        Fbeta(i) = (1+beta^2)*P(i)*R(i)/(beta^2*P(i) + R(i));
    end
end
macrofbeta = mean(Fbeta);

end

