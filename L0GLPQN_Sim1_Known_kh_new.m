function uout = L0GLPQN_Sim1_Known_kh_new(pho, k, h)

addpath(genpath("./PQN"));
%addpath(genpath("/home/huanxi/Documents/Tools/glmnet-matlab"));
%addpath(genpath("/home/huanxi/Documents/Project/L0GL"));
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

load ./Sim1_Known_kh/Xnd.txt
load ./Sim1_Known_kh/yn.txt
load ./Sim1_Known_kh/wd.txt
X = Xnd;
[nInstances, nVars] = size(X);
y = yn;
w = wd;
[rid, cid] = find(wd);
utrue = w;
utrue(rid,1) = 1.;
%k=80;
%h=10;
%pho = 5;%sqrt(nInstances);
gn = 10;
for i =1:gn
    Group{i}=(i-1)*100+1:i*100;
end



% Initial guess of parameters
uSimplex = ones(nVars,1)*(1/nVars);%zeros(nVars,1);

% Set up Objective Function
funObj = @(w)L0Obj(w,X,y,pho);

% Set up Simplex Projection Function
funProj = @(w)ProjCSimplexGL_Gruobi_v2(w,k,Group,h);

% Solve with PQN
options.maxIter = 50;
options.SPGiters = 5;
tStart = cputime;
[uout, err, Timespent, obj] = minConF_PQN_V3(funObj,uSimplex,funProj,options);
tEnd = cputime - tStart

[uout, zout] = ProjCSimplexGL_Gruobi_v2(uout, k, Group,h);
[B, Ranktmp] = sort(-uout);
Rank = sort(Ranktmp(1:k));
uout(Ranktmp(1:k));

Indtrue = find(utrue);
C = intersect(Rank,Indtrue);
length(C)/k


