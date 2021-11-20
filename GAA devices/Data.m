%% TDDB current statistics

load('C:\Users\Mason\Documents\Project\Matlab Project\MMC-GAA\Results\TDDB_2_6_GAA.mat')


figure
hold on
for i=1:length(final_list(:,1))
    s=cell2mat(final_list(i,2));
    u1=s(:,1);
    u2=s(:,3);
    plot(u1,u2,'b')
end

hold off


%% Simulation Clustering
s1=load('C:\Users\Mason\Documents\Project\Matlab Project\MMC-GAA\Results\TDDB_2.1_6_GAA.mat');
s2=load('C:\Users\Mason\Documents\Project\Matlab Project\MMC-GAA\Results\2.1V_8_GAA.mat');
s3=load('C:\Users\Mason\Documents\Project\Matlab Project\MMC-GAA\Results\2.1V_10_GAA.mat');



ss1=cell2mat(s1.final_list(:,1));
data1=sort(ss1(:,1));
ss2=cell2mat(s2.final_list(:,1));
data2=sort(ss2(:,1));
ss3=cell2mat(s3.final_list(:,1));
data3=sort(ss3(:,1));

% ind=round(length(data)/3);
% ff=data(end-ind+1:1:end);
% ff=ff-data(end-ind);
% ff=ff*4/5;
% data(end-ind+1:1:end)=ff+data(end-ind);


n1=1:length(data1);
p1=(n1-0.3)/(length(data1)+0.4);
u1=log(-log(1-p1));

n2=1:length(data2);
p2=(n2-0.3)/(length(data2)+0.4);
u2=log(-log(1-p2));

n3=1:length(data3);
p3=(n3-0.3)/(length(data3)+0.4);
u3=log(-log(1-p3));
% plot(log(data),u,'*')

% [p4,data4,final4,curve4]=Cluster_fitting_bicluster(data,1/3,[],[],[]);
[sss1,data1,final1,curve1]=Cluster_fitting_single(data1,1/3,[],[],[]);
[sss2,data2,final2,curve2]=Cluster_fitting_single(data2,1/3,[],[],[]);
[sss3,data3,final3,curve3]=Cluster_fitting_single(data3,1/3,[],[],[]);

figure
hold on
plot(curve1(1,:),curve1(2,:))
plot(data1,u1,'*')
plot(curve2(1,:),curve2(2,:))
plot(data2,u2,'*')
plot(curve3(1,:),curve3(2,:))
plot(data3,u3,'*')
hold off
%%


load('C:\Users\Mason\Documents\Project\Matlab Project\MMC-GAA\New folder\0.8_162_50_12_PIT2_d8_1.8.mat')
uu=cell(1,length(final_list(:,1)));
filename = '0.8V_IL0.6_HK1.8_D8.xlsx';
factor=1;
figure
hold on
for i=1:length(final_list(:,1))
s=cell2mat(final_list(i,2));
%s(1,1)=8e5;
u1=s(:,1)*factor;
u2=s(:,3);
u2=u2/((u1(2)/factor)^(slope))*(u1(5)^(slope));
u2=u2/3;
ind=find(logical((u2-2e-6)>0),1);
u1(ind+1:end)=[];
u2(ind+1:end)=[];
u2(end)=2e-6;
A=[{'Time','Current'}];
xlswrite(filename,A,i,'A1')
xlswrite(filename,[u1,u2],i,'A2')
plot(u1,u2,'r')
end