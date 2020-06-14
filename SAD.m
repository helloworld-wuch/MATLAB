function sad = SAD(V,newV)
%SAD 此处显示有关此函数的摘要
%   此处显示详细说明
[~,n]=size(V);
s=[];
for i=1:n
    a=V(:,i);
    b=newV(:,i);
    c=acos((a'*b)/(sqrt(a'*a)*sqrt(b'*b)));
    s(:,i)=c;
end
sad=mean(s);
end

