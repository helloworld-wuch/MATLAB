function [W,H,objhistory,iter,elapsed] = NMFSC( V, rdim, sW, sH, fname, showflag, stopconv, tol,timelimit, maxiter )

    
% 检查是否有非负数据
if min(V(:))<0, error('Negative values in data!'); end
    
% 防止溢出
V = V/max(V(:));%缩小到【0，1】

% 数据维度
vdim = size(V,1);
samples = size(V,2);

cons=zeros(samples,samples);  %生成samples*samples的零矩阵
consold=cons;
inc=0;
j=0;
    
% 创建初始矩阵
W = abs(randn(vdim,rdim)); %让基矩阵和丰度矩阵非负
H = abs(randn(rdim,samples));
% H = H./(sqrt(sum(H.^2,2))*ones(1,samples));

%使初始矩阵具有正确的稀疏性
if ~isempty(sW), 
    L1a = sqrt(vdim)-(sqrt(vdim)-1)*sW;  %在欧几里得距离最接近的稀疏约束,计算L1范数（固定L2范数保持基本不变(L2设置为1)，
    %通过调整的L1范数使得向量满足给定稀疏度
    for i=1:rdim, 
        W(:,i) = projfunc(W(:,i),L1a,1,1); %W（：，i）表示W矩阵的第i列
    end
end
if ~isempty(sH), 
    L1s = sqrt(samples)-(sqrt(samples)-1)*sH; 
    for i=1:rdim, 
        H(i,:) = (projfunc(H(i,:)',L1s,1,1))'; 
    end
end

% 初始化显示
if showflag,
    figure(1); clf; % 显示能量和稀疏
    figure(2); clf; % 显示目标函数
    drawnow;        % 刷新屏幕
end

% 计算初始目标
objhistory = 0.5*sum(sum((V-W*H).^2));

% 初始步长
stepsizeW = 0.5;
stepsizeH = 0.5;

timestarted = clock;  %获取当前时间
elapsed = etime(clock,timestarted);%计算函数运行时间
% 开始迭代
iter = 0;
for iter=1:maxiter,
    % 停止条件
    if objhistory(end) < tol | elapsed > timelimit,
        break;
    end

    % Show progress
    fprintf('[%d]: %.5f\n',iter,objhistory(end));    

     % 每隔一段时间保存一次
    if rem(iter,5)==0,
        % 每 5 次迭代测试收敛
        j=j+1;
        
        % 构造连接矩阵
        [y,index]=max(H,[],1);   %查找最大因素
        mat1=repmat(index,samples,1);  % 利差指数下降
        mat2=repmat(index',1,samples); % 点差索引右侧
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
%         fprintf('Saving...');
%         save(fname,'V','W','H','iter','objhistory','elapsed','inc');
%         fprintf('Done!\n');
    end
	
    % Show stats 显示统计信息
    if showflag & (rem(iter,5)==0),
        figure(1);
        subplot(3,1,1); 
        bar(sqrt(sum(W.^2)).*sqrt(sum(H'.^2)));
        
        cursW = (sqrt(vdim)-(sum(abs(W))./sqrt(sum(W.^2))))/(sqrt(vdim)-1);
        subplot(3,1,2); 
        bar(cursW);
        
        cursH = (sqrt(samples)-(sum(abs(H'))./sqrt(sum(H'.^2))))/(sqrt(samples)-1);
        subplot(3,1,3); 
        bar(cursH);
        
        if iter>1,
            figure(2);
            plot(objhistory(2:end));
        end
        drawnow;
    end
    
    %更新迭代计数
    iter = iter+1;    
    
    % Save old values
    Wold = W;
    Hold = H;
        
    % ----- Update H ---------------------------------------
    %只在权重矩阵H上加稀疏约束
    if ~isempty(sH),
        
        % H 的渐变
        dH = W'*(W*H-V);
        begobj = objhistory(end);
        
        % 确保我们降低目标！
        while 1,
            % 向负梯度方向迈出一步，并投影
            Hnew = H - stepsizeH*dH;
            for i=1:rdim, 
                Hnew(i,:) = (projfunc(Hnew(i,:)',L1s,1,1))'; %对H的每一行进行稀疏化
            end

            % 计算新目标
            newobj = 0.5*sum(sum((V-W*Hnew).^2));

            % 如果目标下降，我们可以继续...
            if newobj<=begobj,
                break;
            end

            %...否则减小步长大小，然后重试
            stepsizeH = stepsizeH/2;
            fprintf('.');
            if stepsizeH<1e-200, 
                fprintf('Algorithm converged.\n');
                return; 
            end
        end
        
        % 稍微增加步长
        stepsizeH = stepsizeH*1.2;
        H = Hnew;

    else
        % 使用标准 NMF 多页更新规则进行更新
        H = H.*(W'*V)./(W'*W*H + 1e-9);

        % 重新规范化，使 H 行具有恒定的能量
        norms = sqrt(sum(H'.^2));
        H = H./(norms'*ones(1,samples));
        W = W.*(ones(vdim,1)*norms);  
    end
    
    
    % ----- Update W ---------------------------------------

    if ~isempty(sW),    
        % W 的渐变
        dW = (W*H-V)*H';
        begobj = 0.5*sum(sum((V-W*H).^2));
	
        % 确保我们降低目标！
        while 1,
            % 向负梯度方向迈出一步，并投影
            Wnew = W - stepsizeW*dW;
            norms = sqrt(sum(Wnew.^2));
            for i=1:rdim, 
                Wnew(:,i) = projfunc(Wnew(:,i),L1a*norms(i),(norms(i)^2),1); 
                %对W的每一列进行稀疏化
            end
	
            % 计算新目标
            newobj = 0.5*sum(sum((V-Wnew*H).^2));
	    
            % 如果目标下降，我们可以继续...
            if newobj<=begobj,
                break;
            end
	    
            % ...否则减小步长大小，然后重试
            stepsizeW = stepsizeW/2;
            fprintf(',');
            if stepsizeW<1e-200, 
                fprintf('Algorithm converged.\n');
                return; 
            end
        end
        
        % 稍微增加步长
        stepsizeW = stepsizeW*1.2;
        W = Wnew;
    
    else
        %  使用标准 NMF 多页更新规则进行更新
        W = W.*(V*H')./(W*H*H' + 1e-9);	
    end
    
    % 计算目标
    newobj = 0.5*sum(sum((V-W*H).^2));
    objhistory = [objhistory newobj];
end
