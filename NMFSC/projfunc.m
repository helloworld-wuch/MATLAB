function [v,usediters] = projfunc( s, k1, k2, nn )
%实现 L1\L2 标准投影
%给定向量 s，找到在欧几里距离中最接近 s 的矢量 v 与 sum（abs（v）=k1 和 sum（v^2）_k2。
%如果设置了二进制标志 nn，则向量 v 还被限制为非负 （v=0）。
    
% 问题维度
N = length(s);

% 如果未设置非负性标志，请记录标志并采取 abs
if ~nn,
    isneg = s<0;
    s = abs(s);
end

% 首先将点投影到总约束超平面
v = s + (k1-sum(s))/N;

% 初始化零coeff（最初，没有假定为零元素）
zerocoeff = [];

j = 0;
while 1,

    % This does the proposed projection operator
    midpoint = ones(N,1)*k1/(N-length(zerocoeff));
    midpoint(zerocoeff) = 0;
    w = v-midpoint;
    a = sum(w.^2); %∑（yi-mi）^2*α^2﹢2∑mi（yi-mi）α﹢∑m^2i-l2^2=0 的系数为a,b,c
    b = 2*w'*v;
    c = sum(v.^2)-k2;
    alphap = (-b+real(sqrt(b^2-4*a.*c)))./(2*a);
    %fprintf('alpha1 %.5f ', size(alphap,1));
    v = w*alphap + v;
    
    if all(v>=0),
	%如果v中所有元素均为非负，运算结束。
	usediters = j+1;
	break;
    end
        
    j = j+1;
        
    % 将内数设置为零，从休息中减去适当的数量
    zerocoeff = find(v<=0);
    v(zerocoeff) = 0;
    tempsum = sum(v);
    v = v + (k1-tempsum)/(N-length(zerocoeff));
    v(zerocoeff) = 0;
            
end

% 如果未设置非负性标志，则将符号返回到解决方案
if ~nn,
    v = (-2*isneg + 1).*v;
end

% Check for problems
if max(max(abs(imag(v))))>1e-10,
    error('Somehow got imaginary values!');
end
