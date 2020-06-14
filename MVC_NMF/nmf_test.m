clear all;
clear; clc

% S=load('C:\Users\Archive\Desktop\img\Indian_pines_gt.mat');%
% S1=struct2cell(S);%将结构转变为元胞数组
% V=cell2mat(S1);%将元胞数组转变成为普通的矩阵 
% V= V/max(V(:));
load('Indian_pines_corrected.mat');
X=indian_pines_corrected;
M=X(:,:,140);
V=mapminmax(M,0,1);
mixed=V;
M = 29;
N = 5;

% remove noise
c = test_p(mixed');%通过主成分分析测出端元数
c=int8(c);
[UU, SS, VV] = svds(mixed,c);%SS为c个最大特征值构成的对角矩阵； UU,VV分别是相应的特征值对应的列向量和行向量
%mixed约等于UU*SS*VV'
[L I]=size(mixed);
r_m = mean(mixed,2); %对每行求均值     
R_m = repmat(r_m,[1 I]); % 扩展列向量，变成I个相同列向量构成的矩阵
R_o = mixed - R_m; 
x_p =  UU' * R_o;  
SNR=estimate_snr(mixed,r_m,x_p);%测量信噪比

Lowmixed = UU'*mixed;
mixed = UU*Lowmixed;

% HySime algorithm 
verbose = 'on';
c=double(c);
%[A_vca, EndIdx] = vca(mixed,'Endmembers', c,'verbose','on');
[A_vca, EndIdx] = vca(mixed,'Endmembers', c,'SNR', SNR,'verbose','on');

% 利用全局最小二乘法对丰度矩阵进行演算
warning off;
AA = [1e-5*A_vca;ones(1,length(A_vca(1,:)))];
s_fcls = zeros(length(A_vca(1,:)),M*N);
for j=1:M*N
    r = [1e-5*mixed(:,j); 1];%为了让其必然有解
    % s_fcls(:,j) = nnls(AA,r);
    s_fcls(:,j) = lsqnonneg(AA,r);%求解非负线性最小二乘问题，返回s_fcls(:,j)大于0情况下的最小向量
end

% use vca to initiate
Ainit = A_vca;%基矩阵
sinit = s_fcls;%系数矩阵


% use vca to initiate
PrinComp= pca(mixed');     %返回N的主成分系数 通过P数据矩阵X, X的行对应于观测值，列对应于变量。
meanData = mean(mixed');

tol = 0.00001;
maxiter = 5000;
T = 0.015;%φ=λ/(p-1)!系数
showflag = 1;
figure(2);
imshow(mixed);
[Aest, sest] = mvcnmf_new(mixed,Ainit,sinit,UU,PrinComp,meanData,T,tol,maxiter,showflag,2,1);
newAS=Aest*sest;
figure(3);
imshow(newAS);
SAD_4=sum(sum(abs(V-newAS)))/(145*145);
RMSE_4=sum(sqrt(mean((V-newAS).^2)))/145;
subplot(1,2,1);imshow(V);title('Original');
subplot(1,2,2);imshow(newAS);title('MVCNMF');xlabel("SAD:"+SAD_4+" RMSE:"+RMSE_4);
