clear;
clc;
load('Indian_pines_corrected.mat');
indian=indian_pines_corrected;
addpath(genpath('C:\Users\Archive\Documents\MATLAB'));
X=indian(:,:,140);
V=mapminmax(X,0,1);
figure(1);
axis([0 200 0 0.4]);
drawnow;
X_rand=[];
j=0;
for i=2:10:200
    j=j+1;
    X=indian(:,:,i);
    V=mapminmax(X,0,1);
    X_rand(1,j)=i;
    
    %VCA
    [newV_1]=plot_VCA(V);
    SAD_1=SAD(V,newV_1);
    RMSE_1=sum(sqrt(mean((V-newV_1).^2)))/145;
    Y_SAD1(1,j)=SAD_1;
    Y_RMSE1(1,j)=RMSE_1;
    
    %NMF
    [X_2]=NMF(V);
    Y_SAD2(1,j)=SAD(V,X_2);
    Y_RMSE2(1,j)=sum(sqrt(mean((V-X_2).^2)))/145;
    
    %NMFSC
    [A_3,S_3]=plot_NMFSC(V);
    newV_3=A_3*S_3;
    RMSE_3=sum(sqrt(mean((V-newV_3).^2)))/145;
    Y_SAD3(1,j)=SAD(V,newV_3);
    Y_RMSE3(1,j)=RMSE_3;
    
    %MVCNMF
    [A_4,S_4]=plot_MVCNMF(V);
    newV_4=A_4*S_4;
    SAD_4=SAD(V,newV_4);
    RMSE_4=sum(sqrt(mean((V-newV_4).^2)))/145;
    Y_SAD4(1,j)=SAD_4;
    Y_RMSE4(1,j)=RMSE_4;
    
    %WECNMF
    [A_5,S_5]=plot_WECNMF(V);
    newV_5=A_5*S_5;
    SAD_5=SAD(V,newV_5);
    RMSE_5=sum(sqrt(mean((V-newV_5).^2)))/145;
    Y_SAD5(1,j)=SAD_5;
    Y_RMSE5(1,j)=RMSE_5;
    
    %GNMF
    [newV_6]=plot_GNMF(V);
%     newV_6=A_6*S_6;
    SAD_6=SAD(V,newV_6);
    RMSE_6=sum(sqrt(mean((V-newV_6).^2)))/145;
    Y_SAD6(1,j)=SAD_6;
    Y_RMSE6(1,j)=RMSE_6;
end
% plot(X_rand,Y_SAD1,'r',X_rand,Y_SAD2,'g',X_rand,Y_SAD3,'b',X_rand,Y_SAD4,'y',X_rand,Y_SAD5,'k',X_rand,Y_SAD6,'m');
% plot(X_rand,Y_RMSE1,'r',X_rand,Y_RMSE2,'g',X_rand,Y_RMSE3,'b',X_rand,Y_RMSE4,'y',X_rand,Y_RMSE5,'k',X_rand,Y_RMSE6,'m');
% xlabel("rand");
% ylabel("SAD");% ylabel("RMSE");
% legend('VCA','NMF','NMFSC','MVCNMF','WECNMF','GNMF');
% title('SAD比较');
% title('RMSE比较');