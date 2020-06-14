function [A,S] = plot_NMFSC(V)
rdim = 100;
sW = 0.6;%初始稀疏度
sH = 0.6;
fname = ['Indian_pines_NMFSC'];
showflag = 0;%显示图像
tol = 0.00001;%精度
stopconv = 30;%停止状态
timelimit = 100;%时间限制，太长了也要停止，但结果无效
maxiter = 5000;%迭代次数
[A,S] = NMFSC( V, rdim, sW, sH, fname, showflag, stopconv, tol,timelimit, maxiter );
end

