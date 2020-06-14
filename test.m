clear;
clc;
% S=load('C:\Users\Archive\Desktop\img\Indian_pines_gt.mat');%
load('Indian_pines_corrected.mat');
indian=indian_pines_corrected;
% V=reshape(V,4205,1000);
X=indian(:,:,140);
V=mapminmax(X,0,1);
% S1=struct2cell(S);%将结构转变为元胞数组
% V=cell2mat(S1);%将元胞数组转变成为普通的矩阵 
%V = V/max(V(:));
addpath(genpath('C:\Users\Archive\Documents\MATLAB'));
value=input('输入选择 1.VCA  2.NMF   3.NMFSC  4.MVCNMF   5.WECNMF  6.GNMF>>>>');
switch(value)
    case 1
        [newV_1]=plot_VCA(V);
        subplot(3,3,2);imshow(newV_1);title('VCA');
    case 2
        [X_2]=NMF(V);
        subplot(3,3,3);imshow(X_2);title('NMF');
    case 3
        [A_3,S_3]=plot_NMFSC(V);
        newV_3=A_3*S_3;
        subplot(3,3,4);imshow(newV_3);title('NMFSC');
    case 4
        [A_4,S_4]=plot_MVCNMF(V);
        newV_4=A_4*S_4;
        subplot(3,3,5);imshow(newV_4);title('MVCNMF');
    case 5
        [A_5,S_5]=plot_WECNMF(V);
        newV_5=A_5*S_5;
        subplot(3,3,6);imshow(newV_5);title('WECNMF');
    case 6
        [newV_6]=plot_GNMF(V);
        subplot(3,3,7);imshow(newV_6);title('GNMF');
end
