f=load('Reuters.mat');%=load('C:\Users\User\Desktop\research\multiviewlearning\shiguoxin\bbc_seg14of4.mat');
data=f.data;
label=f.labels;
% addpath('C:\Users\User\Desktop\research\kernelclusteringexp')
para1=[.01 .1 1];
para2=[100,500,1000,1500,2000];
para3=[.01];

for i=1:size(data,1)
dist = max(max(data{i})) - min(min(data{i}));
m01 = (data{i} - min(min(data{i})))/dist;
data{i} = 2 * m01 - 1;
end

for i=1:length(para1)
    for j=1:length(para2)
        for k=1:length(para3)
            result=multigraph(data,label,para1(i),para2(j),para3(k))
            dlmwrite('reuters.txt',[para1(i) para2(j) para3(k) result(1,:) result(2,:) result(3,:)   ],'-append','delimiter','\t','newline','pc');
        end
    end
end
        
