function [X, SAD, RMSE] = NMF(V)
V = V/max(V(:));
[n1 n2]=size(V);
r=100;
W=rand(n1,r);
H=rand(r,n2);
maviter=500;
for iter=1:maviter
   W=W.*(V*H')./(W*H*H'+eps);
   H=H.*(W'*V)./((W')*W*H+eps);
end
X=W*H;  
return;

