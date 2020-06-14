function [pi,P]= Find_P(S)
%权矩阵P即丰度矩阵行和的倒数
[x,y]=size(S);
pi=sum(S,2);
P=zeros(x,y);
for i=1:x
    P(i,i)=1/pi(i,1);  
end
return

