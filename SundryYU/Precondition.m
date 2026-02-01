function data = Precondition(X,style)
% Input:
%       - X: n*d的数据
%       - style: 数据预处理的方式;
%           - 0: 标准化
%           - 1: 将各特征压缩至单位球
%           - 2: 将各样本压缩至单位球

% Output:
%       - data: n*d的，预处理之后的数据

if nargin <= 1
    style = 0;
end

% 对样本进行标准化处理
if style == 0 % 对数据进行标准化
    % 计算每个特征的均值和标准差
    mu = mean(X);
    sigma = std(X);
    % 对每个样本进行标准化处理
    data = (X - mu) ./ sigma;
elseif style == 1 % 将每个特征压缩至单位球
    l2 = X.^2;
    sl2 = sum(l2,style).^(0.5);
    data = X./sl2;
else % 将每个样本压缩至单位球
    l2 = X.^2;
    sl2 = sum(l2,style).^(0.5);
    data = X./sl2;
end

end

