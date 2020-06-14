clear all;
clc;
load('Indian_pines_corrected.mat');
X=indian_pines_corrected;
M=X(:,:,160);
V=mapminmax(M,0,1);
%imshow(M);
[n1 n2]=size(V);
r=120;
W=rand(n1,r);
H=rand(r,n2);
maviter=40;
for iter=1:maviter
   W=W.*(V*H')./(W*H*H'+eps);
   H=H.*(W'*V)./((W')*W*H+eps);
end
d=W*H;  
figure(2);
imshow(d);
SAD=sum(acos(sum(V'*d)/(sqrt((V.^2)*(d).^2))))/145;%光谱角距离
RMSE=sum(sqrt(mean((V-d).^2)))/145;%均方根误差
