clear ;
load('Indian_pines_corrected.mat');
X=indian_pines_corrected;
M=X(:,:,160);
M=mapminmax(M,0,1);
imshow(M);