function [A,S]=wecnmf(X,k)

maxiter=1000;
[r,c]=size(X); 
A=rand(r,k);
S=rand(k,c);
showflag=0;
objhistory = 0.5*sum(sum((X-A*S).^2));%目标函数
timestarted = clock;  %获取当前时间
elapsed = etime(clock,timestarted);%计算函数运行时间
timelimit=100;
tol=0.0001;
cons=zeros(c,c);  %生成samples*samples的零矩阵
consold=cons;
inc=0;
stopconv=42;
fname=['Indian_pines'];
% 初始化显示
if showflag,
    figure(1); clf; % 显示能量和稀疏
    figure(2); clf; % 显示目标函数
    drawnow;        % 刷新屏幕
end


for iter=1:maxiter
     % 停止条件
    if objhistory(end) < tol || elapsed > timelimit
        break;
    end
    
    fprintf('[%d]: %.5f\n',iter,objhistory(end));    
    
     % 每隔一段时间保存一次
    if rem(iter,5)==0
        % 每 5 次迭代测试收敛   
        % 构造连接矩阵
        [~,index]=max(S,[],1);   %查找最大因素
        mat1=repmat(index,c,1);  % 利差指数下降
        mat2=repmat(index',1,c); % 点差索引右侧
        cons=mat1==mat2;
        
        if(sum(sum(cons~=consold))==0) % 连接矩阵未更改
            inc=inc+1;                     %累积计数
        end
        fprintf('\t%d\t%d\t%d\n',iter,inc,sum(sum(cons~=consold))),
        
        if(inc>=stopconv)
            break,                % 假设融合是连接停止变化 
        end
        
        consold=cons;
        
        elapsed = etime(clock,timestarted);
        fprintf('Saving...');
        save(fname,'X','A','S','iter','objhistory','elapsed','inc');
        fprintf('Done!\n');
    end
    
    %显示迭代信息
    if showflag & (rem(iter,5)==0),
        figure(1);
        subplot(3,1,1); 
        bar(sqrt(sum(A.^2)).*sqrt(sum(S'.^2)));

        cursA = (sqrt(r)-(sum(abs(A))./sqrt(sum(A.^2))))/(sqrt(r)-1);
        subplot(3,1,2); 
        bar(cursA);

        cursS = (sqrt(c)-(sum(abs(S'))./sqrt(sum(S'.^2))))/(sqrt(c)-1);
        subplot(3,1,3); 
        bar(cursS);

    if iter>1,
        figure(2);
        plot(objhistory(2:end));
    end
    drawnow;
    end
    
     
    %迭代内容
    [pi,P]=Find_P(S);%权矩阵P即丰度矩阵行和的倒数
    M_m=mean(X,2); %求各列向量的均值
    M= repmat(M_m,[1 k]); % 扩展到1到N行
    WED=0.5*sum(sum(((A-M)*P).^2));
    
    WED_A=(A-M)*P*(P'); %对基矩阵求偏导
    for j=1:k %对系数矩阵求偏导
        AM=A(:,j)-M(:,j);
        WED_S(j,:)=-(1/power(pi(j,1),3))*power(AM,2);
    end

    A=A.*(X*S'-0.5*WED_A)./((A*S*(S'))+eps);
    S=S.*((A')*X-0.5*WED_S)./(((A')*A*S)+eps);
    
    % 计算目标
    newobj = 0.5*sum(sum((X-A*S).^2))+WED;
    objhistory = [objhistory newobj];
end

return
