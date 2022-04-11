uout = L0GLPQN_IMATNIB(0.2, 50, 10)

load ./RealWorld/IMATINIB_X_Train.txt
load ./RealWorld/IMATINIB_y_Train.txt
X = IMATINIB_X_Train;
[nInstances, nVars] = size(X);
y = IMATINIB_y_Train;

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
gn = length(Group);
B = zeros(nVars, gn);
for i =1:gn
    B(Group{i},i)=1;
end

[val, rank] = sort(gout,'descend');

Indselect = zeros(nVars,1);
for i=1:10
    Indselect(Group{rank(i)}) = 1;
end
