function uout = L0GLPQN_IMATNIB(pho, k, h)

addpath(genpath("/home/huanxi/Documents/Project/PQN/PQN"));
addpath(genpath("/home/huanxi/Documents/Tools/glmnet-matlab"));
addpath(genpath("/home/huanxi/Documents/Project/L0GL"));
addpath(genpath("/home/huanxi/Documents/Tools/gurobi910/linux64/matlab/"));


% Generate Syntehtic Data
% nInstances = 200:100:1200;
% nVars = 5000;
% k = 20;
% pho=0.2;
% SNR = 6;
% 
% T = toeplitz(pho.^(1:nVars));
% RandIndtmp = randperm(nVars);
% RandInd = sort(RandIndtmp(1:k));
% RandInt = 2*randi([0,1],k,1)-1;
% X =  mvnrnd(zeros(1,nVars), T, nInstances(i));
% w = zeros(nVars,1);
% w(RandInd) = RandInt;%rand(nVars,1).*(rand(nVars,1) > .5);
% utrue = w;% Set up Simplex Projection Function
% 
% utrue(RandInd) = 1.0;
% NoiseT = normrnd(0,1,nInstances(i),1);
% Noise = ((norm(X*w)/sqrt(SNR))/norm(NoiseT))*NoiseT;
% y = X*w + Noise;
% pho = sqrt(nInstances(i));

load ./RealWorld/IMATINIB_X_Train.txt
load ./RealWorld/IMATINIB_y_Train.txt
X = IMATINIB_X_Train;
[nInstances, nVars] = size(X);
y = IMATINIB_y_Train;
% [rid, cid] = find(wd);
% utrue = w;
% utrue(rid,1) = 1.;
% %k=80;
%h=10;
%pho = 5;%sqrt(nInstances);
%gn = 50;
% for i =1:gn
%     Group{i}=(i-1)*10+1:i*10;
% end
utrue = rand(nVars,1);

fid = fopen('./RealWorld/Pathway_Group_pathway7.txt');
line1 = fgetl(fid);
res=line1;
while ischar(line1)
line1 = fgetl(fid);
res =char(res,line1);
end
fclose(fid);
for i=1:size(res,1)-1
  Group{i}=str2num(res(i,:))+1;
end



% Initial guess of parameters
uSimplex = ones(nVars,1)*(1/nVars);%zeros(nVars,1);

% Set up Objective Function
funObj = @(w)L0Obj(w,X,y,pho);

% Set up Simplex Projection Function
funProj = @(w)ProjCSimplexGL_Gruobi(w,k,Group,h);

% Solve with PQN
options.maxIter = 50;
options.SPGiters = 3;
tStart = cputime;
[uout, err, Timespent, obj] = minConF_PQN_V2(funObj,uSimplex,funProj,options);
tEnd = cputime - tStart

[uout, zout] = ProjCSimplexGL_Gruobi_v2(uout, k, Group,h);
[B, Ranktmp] = sort(-uout);
Rank = sort(Ranktmp(1:k));
uout(Ranktmp(1:k));

Indtrue = find(utrue);
C = intersect(Rank,Indtrue);
length(C)/k


