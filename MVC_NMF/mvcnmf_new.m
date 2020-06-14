function [A,S,time] = mvcnmf_new(X,Ainit,Sinit,UU,PrinComp,meanData,T,tol,maxiter,showflag,type_alg_S,type_alg_A)

% A,S: output solution
% Ainit,Sinit: initial solutions
% Atrue: true endmembers
% UU: principle components for visualization (SVD)
% PrinComp: principal components for calculating volme (PCA)
% meanData: for calculating volume
% T: annealing temprature
% tol: tolerance for a relative stopping condition
% maxiter: limit of iterations
% showflag: display scatter plot (1)
% type_alg_S: algorithms for estimating S
% type_alg_A: algorithms for estimating A

A = Ainit; S = Sinit; 

% dimensions
c = size(S,1);     % number of endmembers
N = size(S,2);     % number of pixels

% precalculation for visualization
%EM = UU'*Atrue;         % low dimensional endmembers
LowX = UU'*X;           % low dimensional data
inc=0;
j=0;

% PCA to calculate the volume of true EM
% E = [ones(1,c);PrinComp(:,1:c-1)'*(Atrue-meanData'*ones(1,c))];
% vol_t = 1/factorial(c-1)*abs(det(E)); % the volume of true endmembers
% vol = [];

% calculate volume of estimated A
C = [ones(1,c); zeros(c-1,c)];
B = [zeros(1,c-1); eye(c-1)];
Z = C+B*(PrinComp(:,1:c-1)'*(A-meanData'*ones(1,c))); 
detz2 = det(Z)*det(Z);

% one time draw
if showflag,
    startA = UU'*Ainit;
    figure(1),
    for i=1:3
        for j=i+1:3
            subplot(2,2,(i-1)*2+j-i),
            plot(LowX(i,1:6:end),LowX(j,1:6:end),'rx'); %每6给为间隔提取第i行的数据，并用x绘制
            hold on, 
            plot(startA(i,:),startA(j,:),'bo'); %original estimateion
%             plot(EM(i,:), EM(j,:),'go','markerfacecolor','g'); %true endmember
        end
    end
    
end


%


% 计算初始梯度
gradA = A*(S*S') - X*S' + T*detz2*PrinComp(:,1:c-1)*B'*pinv(Z)'; 
gradS = (A'*A)*S - A'*X;
initgrad = norm([gradA; gradS'],'fro'); 
fprintf('Init gradient norm %f\n', initgrad); 
tolA = max(0.001,tol)*initgrad; tolS = tolA;

% 计算初始目标
objhistory = 0.5*sum(sum((X-A*S).^2));
Ahistory = [];
% 初始化显示
if showflag
    figure(1); clf; % 显示能量和稀疏
    figure(2); clf; % 显示目标函数
    drawnow;        % 刷新屏幕
end
% 初始步长
stepsizeW = 1;
stepsizeH = 1;

timestarted = clock;  %获取当前时间
elapsed = etime(clock,timestarted);%计算函数运行时间
% count the number of sucessive increase of obj

timelimit=100;
vdim = size(X,1);
samples = size(X,2);
cons=zeros(samples,samples);  %生成samples*samples的零矩阵
consold=cons;
stopconv=42;
fname=['C:\Users\Archive\Desktop\img\Indian_pines_gt.mat.mat'];
for iter=1:maxiter
    % 停止条件
    if objhistory(end) < tol || elapsed > timelimit
        break;
    end

    % Show progress
    fprintf('[%d]: %.5f\n',iter,objhistory(end));    

     % 每隔一段时间保存一次
    if rem(iter,5)==0
        % 每 5 次迭代测试收敛
        j=j+1;
        
        % 构造连接矩阵
        [~,index]=max(S,[],1);   %查找最大因素
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
        fprintf('Saving...');
%         save(fname,'X','A','S','iter','objhistory','elapsed','inc');
        fprintf('Done!\n');
    end
    
    %更新迭代计数
    iter = iter+1;    
    
    % Save old values
%     Wold = W;
%     Hold = H;
        
    E = [ones(1,c);PrinComp(:,1:c-1)'*(A-meanData'*ones(1,c))];
    vol_e = 1/factorial(c-1)*abs(det(E));
    fprintf('[%d]: %.5f\t',iter,objhistory(end));    
    fprintf('Temperature: %f \t', T);
%     fprintf('Actual Vol.: %f \t Estimated Vol.: %f\n', vol_t, vol_e);
    vol(iter+1) = vol_e;
    
    % real time draw
    if showflag,
        est = UU'*A;      
        Ahistory = [Ahistory est];
        figure(1),
        for i=1:3
            for j=i+1:3
                subplot(2,2,(i-1)*2+j-i),
                plot(est(i,:),est(j,:),'yo'); %estimation from nmf
            end
        end
        drawnow;
    end

    % to consider the sum-to-one constraint
    tX = [X; 20*ones(1,N)];
    tA = [A; 20*ones(1,c)];
        
    % 找到S
    switch type_alg_S
        
        case 1 % 共轭梯度学习
            
            no_iter = 50; 
            S = conjugate(X,A,S,no_iter,PrinComp(:,1:c-1),meanData,T);
            
        case 2 % 梯度下降
            
            tolS = 0.0001;
            [S,gradS,iterS] = steepdescent(tX,tA,S,tolS,200,PrinComp(:,1:c-1),meanData,T,'S');
            if iterS==1,
                tolS = 0.1 * tolS; 
            end
    end
    
        % 找到A   
    switch type_alg_A
        
        case 1 % conjugate gradient learning
            
            no_iter = 50; 
            A = conjugate(X',S',A',no_iter,PrinComp(:,1:c-1),meanData,T);
            A = A';
             
        case 2 % steepest descent
                
            tolA = 0.0001;
            [A,gradA,iterA] = steepdescent(X',S',A',tolA,100,PrinComp(:,1:c-1),meanData,T,'A'); 
            A = A'; gradA = gradA';
            if iterA==1,
                tolA = 0.1 * tolA;
            end
             
    end 
    
    % 计算目标
    newobj = 0.5*sum(sum((X-A*S).^2));
    objhistory = [objhistory newobj];
end

    