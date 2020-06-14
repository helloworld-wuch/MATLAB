function [S,grad,iter] = steepdescent(X,A,S,tol,maxiter,U,meanData,tao, mode)

% S, grad:输出解和梯度
% iter: #iterations used
% X, A: 保持不变，找到S
% tol: 停止的容忍度
% maxiter: limit of iterations
% U, meanData: 主成分和均值数据来计算体积
% tao: 正则化参数

[L,N] = size(X); [c,N] = size(S);

%仅对A约束
cons = 0;
% if L>N
if (mode=='A') % LiHe    
    cons = 1;
    % 体积约束的预先计算
    meanData = meanData'*ones(1,c);
    C = [ones(1,c); zeros(c-1,c)];
    B = [zeros(1,c-1); eye(c-1)];
    
end

% 预计算，减少计算复杂度
AtX = A'*X;
AtA = A'*A; 

alpha = 1; beta = 0.1; sigma = 0.01;
% S = max( temp_T - lambda/mu,0) + min( temp_T + lambda/mu,0); 

for iter=1:maxiter,  
    
    % constraint on S^T
    if cons == 1
        Z = C+B*U'*(S'-meanData);
        ZD = pinv(Z)*B*U';
        detz2 = det(Z)^2;%Z的行列式值（|Z|）的平方
        f = sum(sum((X-A*S).^2)) + tao*det(Z)^2;%目标函数
    end
    
   
    if cons == 1 % 对A的梯度
        grad = AtA*S - AtX + tao*detz2*ZD;
    else
%         grad = AtA*S - AtX + 0.0001*sign(S);
        
        %modified
        grad = AtA*S - AtX;  % 对S的梯度
    end
    
    projgrad = norm(grad(grad < 0 | S >0));
    if projgrad < tol,
        break
    end
           
    % 搜索步长
    for inner_iter=1:50,
        Sn = max(S - alpha*grad, 0); d = Sn-S; %Sn的max将负数赋值为0
        
        if cons == 1
            fn = sum(sum((X-A*Sn).^2)) + tao*det(C+B*U'*(Sn'-meanData))^2;
            suff_decr = fn - f <= sigma*sum(sum(grad.*d));
        else       
            gradd=sum(sum(grad.*d)); dQd = sum(sum((AtA*d).*d));
            suff_decr = 0.99*gradd + 0.5*dQd < 0;
        end
        
        if inner_iter==1, % the first iteration determines whether we should increase or decrease alpha
            decr_alpha = ~suff_decr; Sp = S;
        end
        if decr_alpha, 
            if suff_decr,
                S = Sn; break;
            else
                alpha = alpha * beta;
            end
        else
            if ~suff_decr | Sp == Sn,
                S = Sp; break;
            else
                alpha = alpha/beta; Sp = Sn;
            end
        end
    end
end