function [Ae, Rp, SAD, RMSE] = VCA(R,varargin)
%R为高光谱数据集转换的矩阵 varargin为可变参数列表 
verbose = 'on';                          % 显示有关如何隐蔽警告的消息。
snr_input = 0;                           % default this flag is zero,噪声比
% which means we estimate the SNR
dim_in_par = length(varargin);           %求可变参数个数
if (nargin - dim_in_par)~=1              %nargin为输入变量个数
    error('Wrong parameters');
elseif rem(dim_in_par,2) == 1
    error('Optional parameters should always go by pairs');
else
    for i = 1 : 2 : (dim_in_par-1)       %以1为初始 每次加2 截至到-1 不存在，所以不执行
        switch lower(varargin{i})
          case 'verbose'
               verbose = varargin{i+1};
          case 'endmembers'     
               p = varargin{i+1};
          case 'snr'     
               SNR = varargin{i+1};
               snr_input = 1;       % flag meaning that user gives SNR       
          otherwise
               fprintf(1,'Unrecognized parameter:%s\n', varargin{i});
        end %switch
    end %for
end %if
if isempty(R)
    error('there is no data');
else
    [L, N]=size(R);  % L number of bands (channels)  像元或像素点
                % N number of pixels (LxC)           线
end

p=test_p(R);
p=ceil(p);%向上取整
if (p<0 || p>L || rem(p,1)~=0)
    error('ENDMEMBER parameter must be integer between 1 and L');
end

if snr_input==0
    r_m = mean(R,2);         %求各行向量的均值
    R_m = repmat(r_m,[1 N]); % 扩展到1到N列
    R_o = R - R_m;           % data with zero-mean 
    [Ud,Sd,Vd] = svds(R_o*R_o'/N,p);  %  返回 p 个最大奇异值。
    x_p =  Ud' * R_o;                 % 将零均值数据投影到 p 子空间上

    SNR = estimate_snr(R,r_m,x_p);    %计算图像信噪比

    if strcmp (verbose, 'on')
        fprintf(1,'SNR estimated = %g[dB]\n',SNR);
    end
else   
    if strcmp (verbose, 'on')
        fprintf(1,'input    SNR = %g[dB]\t',SNR); 
    end
end

SNR_th = 15 + 10*log10(p);

if SNR < SNR_th    %PCA降维
    if strcmp (verbose, 'on')
        fprintf(1,'... Select the projective proj.\n',SNR); 
    end

    d = p-1;
    if snr_input==0
         Ud= Ud(:,1:d);      % 将零均值数据投影到 p-1 子空间上
    else
         r_m = mean(R,2);      
         R_m = repmat(r_m,[1 N]); % mean of each band
         R_o = R - R_m;           % data with zero-mean 

         [Ud,Sd,Vd] = svds(R_o*R_o'/N,d);  % computes the p-projection matrix 

         x_p =  Ud' * R_o;                 % project thezeros mean data onto p-subspace

    end

    Rp =  Ud * x_p(1:d,:) + repmat(r_m,[1 N]);      % again in dimension L

    x = x_p(1:d,:);             %  x_p =  Ud' * R_o; is on a p-dim subspace
    c = max(sum(x.^2,1))^0.5;
    y = [x ; c*ones(1,N)] ;     %
else   %SVD降维
    if strcmp (verbose, 'on')
        fprintf(1,'... Select proj. to p-1\n',SNR); 
    end

    d = p;
    [Ud,Sd,Vd] = svds(R*R'/N,d);         % computes the p-projection matrix 

    x_p = Ud'*R;
    Rp =  Ud * x_p(1:d,:);      % again in dimension L (note that x_p has no null mean)

    x =  Ud' * R;
    u = mean(x,2);        %equivalent to  u = Ud' * r_m
    y =  x./ repmat( sum( x .* repmat(u,[1 N]) ) ,[d 1]);  %将x投影到超平面上得到y

end
indice = zeros(1,p);
A = zeros(p,p);
A(p,1) = 1;

for i=1:p
      w = rand(p,1);             %随机生成一个零均值高斯随机向量w；
      f = w - A*pinv(A)*w;
      f = f / sqrt(sum(f.^2));   %正交与子空间的向量；
      
      v = f'*y;                  %将数据y投影到f向量上
      [v_max, indice(i)] = max(abs(v));     %求该投影的极值对应的像元位置
      A(:,i) = y(:,indice(i));              %将结果求取一个投影方向
end
Ae = Rp(:,indice);               %端元的光谱曲线

% SAD=sum(sum(abs(R-Rp)))/(145*145);%光谱角距离
% RMSE=sum(sqrt(mean((R-Rp).^2)))/145;%均方根误差

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of the vca function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function snr_est = estimate_snr(R,r_m,x)  %测试噪声比

         [L , ~]=size(R);           % L number of bands (channels)
                                  % N number of pixels (Lines x Columns) 
         [p, N]=size(x);           % p number of endmembers (reduced dimension)

         P_y = sum(R(:).^2)/N;
         P_x = sum(x(:).^2)/N + r_m'*r_m;
         snr_est = 10*log10( (P_x - p/L*P_y)/(P_y- P_x) );
return;

