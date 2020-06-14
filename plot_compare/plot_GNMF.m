function [newX] = plot_GNMF(X)
%PLOT_GNMF 此处显示有关此函数的摘要
%   此处显示详细说明
nClass = 100;

options = [];
options.WeightMode = 'Binary'; 
W = constructW(X,options);%构造权矩阵，度量向量间贴近程度
options.maxIter = 1000;
options.nRepeat = 1;
options.alpha = 1;
%rand('twister',10);
[U,V,n,e] = GNMF(X',nClass,W,options);
newX=V*U';
end

