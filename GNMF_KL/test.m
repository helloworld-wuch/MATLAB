clear;
load('TDT2_all.mat');
idx = gnd <= 10;
fea = fea(idx,:);
gnd = gnd(idx);
nClass = length(unique(gnd));

%tfidf weighting and normalization 
fea = tfidf(fea);

disp('GNMF_KL...');	
Woptions = [];
Woptions.WeightMode = 'Cosine';
Woptions.k = 7;
W = constructW(fea,Woptions);

GNMFKLoptions = [];
GNMFKLoptions.maxIter = 50;
GNMFKLoptions.alpha = 100;
GNMFKLoptions.weight = 'NCW';
nFactor = 10;
rand('twister',5489);
[U, V] = GNMF_KL(fea', nFactor, W, GNMFKLoptions); %'

rand('twister',5489);
label = litekmeans(V,nClass,'Replicates',10); %'
label = bestMap(gnd,label);
AC = length(find(gnd == label))/length(gnd);
MIhat = MutualInfo(gnd,label);
disp(['Clustering in the GNMF_KL space. AC=',num2str(AC),' NMI=',num2str(MIhat)]);