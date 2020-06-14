clear;
clc;
% S=load('C:\Users\Archive\Desktop\img\Indian_pines_gt.mat');%
% S1=struct2cell(S);%将结构转变为元胞数组
% V=cell2mat(S1);%将元胞数组转变成为普通的矩阵 
% V = V/max(V(:));
load('Indian_pines_corrected.mat');
X=indian_pines_corrected;
M=X(:,:,160);
V=mapminmax(M,0,1);
k=100;
[A,S]=wecnmf(V,k);
X=A*S;
% imshow(X);
SAD_5=sum(sum(abs(V-X)))/(145*145);
RMSE_5=sum(sqrt(mean((V-X).^2)))/145;