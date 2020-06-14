clear all;
clc;
load('Indian_pines_corrected.mat');
indian=indian_pines_corrected;
X=indian(:,:,140);
V=mapminmax(X,0,1);%归一化
rdim = 100;
sW = 0.6;%初始稀疏度
sH = 0.6;
fname = ['Indian_pines'];
showflag = 0;%显示图像
tol = 0.00001;%精度
stopconv = 30;%停止状态
timelimit = 100;%时间限制，太长了也要停止，但结果无效
maxiter = 5000;%迭代次数
[W,H,objhistory,iter,elapsed] = NMFSC( V, rdim, sW, sH, fname, showflag, stopconv, tol,timelimit, maxiter );
WH=W*H;
% V = V/max(V(:));
SAD=SAD(V,WH);%光谱角距离

RMSE=sum(sqrt(mean((V-WH).^2)))/145;