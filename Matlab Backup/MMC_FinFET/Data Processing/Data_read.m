[Type,Sheet,Format]=xlsfinfo('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05102016\TDDB_125_fr272829\125C\CVS#2@2.xls');
Sheet(2:3)=[];
res={};
figure;
hold on
for i=1:length(Sheet)
    s=xlsread('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05102016\TDDB_125_fr272829\125C\CVS#2@2.xls',Sheet{i},'A:B');
    s1=zeros(length(s(:,1))+1,2);
    s1(1,1)=0.3;
    s1(1,2)=s(1,2);
    s1(2:end,:)=s;
    res(1,i)={s1};
    %plot(res{1,i}(:,1),res{1,i}(:,2),'g');
end


[Type2,Sheet2,Format2]=xlsfinfo('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05102016\TDDB_125_fr272829\125C\CVS#2.1@2.1.xls');
Sheet2(2:3)=[];
for j=1:length(Sheet2)
    ss=xlsread('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05102016\TDDB_125_fr272829\125C\CVS#2.1@2.1.xls',Sheet2{j},'A:B');
    s2=zeros(length(ss(:,1))+1,2);
    s2(1,1)=0.3;
    s2(1,2)=ss(1,2);
    s2(2:end,:)=ss;
    res(2,j)={s2};
    %plot(res{2,j}(:,1),res{2,j}(:,2),'r');
end

[Type3,Sheet3,Format3]=xlsfinfo('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05102016\TDDB_125_fr272829\125C\CVS#2.2@2.2.xls');
Sheet3(2:3)=[];
for k=1:length(Sheet3)
    sss=xlsread('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05102016\TDDB_125_fr272829\125C\CVS#2.2@2.2.xls',Sheet3{k},'A:B');
    s3=zeros(length(sss(:,1))+1,2);
    s3(1,1)=0.3;
    s3(1,2)=sss(1,2);
    s3(2:end,:)=sss;
    res(3,k)={s3};
    %plot(res{3,k}(:,1),res{3,k}(:,2),'b');
end


[Type4,Sheet4,Format4]=xlsfinfo('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05192016\CVS#FR7-9_125C_2.05V@1.xls');
Sheet4(2:3)=[];
for k=1:length(Sheet4)
    sss_4=xlsread('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05192016\CVS#FR7-9_125C_2.05V@1.xls',Sheet4{k},'A:B');
    s4=zeros(length(sss_4(:,1))+1,2);
    s4(1,1)=0.3;
    s4(1,2)=sss_4(1,2);
    s4(2:end,:)=sss_4;
    res(4,k)={s4};
    %plot(res{4,k}(:,1),res{4,k}(:,2),'y');
end

[Type5,Sheet5,Format5]=xlsfinfo('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05192016\CVS#FR7-9_125C_2.15V@1.xls');
Sheet5(2:3)=[];
for k=1:length(Sheet5)
    sss_5=xlsread('C:\E\Work Station\My Work\Data\IMEC\FinFET\FinFET 05192016\CVS#FR7-9_125C_2.15V@1.xls',Sheet5{k},'A:B');
    s5=zeros(length(sss_5(:,1))+1,2);
    s5(1,1)=0.3;
    s5(1,2)=sss_5(1,2);
    s5(2:end,:)=sss_5;
    res(5,k)={s5};
    %plot(res{5,k}(:,1),res{5,k}(:,2),'m');
end

hold off

list={};
list1={};
result={};
for i=1:5
    for j=1:length(res(i,:))
        if isempty(res{i,j})
            continue
        end
        test=diff(log([res{i,j}]),1);
        test1=diff(test,1);
        o=test(:,2)./test(:,1);
        o1=test1(:,2)./test1(:,1);
        s=o(o>=0);
        s1=o1(o1>=0);
        list(i,j)={s};
        list1(i,j)={s1};
        result(i,j)={[mean([list{i,j}]),var([list{i,j}]),mean([list1{i,j}]),var([list1{i,j}])]};
    end
end

%% Fourier analysis
Data=load('C:\D\Program Files\Project\Matlab Project\Clustering data processing\Fourier\current.mat');
D=Data.res;
list=cell(size(D));
spectrum_list=cell(size(D));

base_f=1;
base_p=[1e-12,1e-11,1e-10,5e-12,3e-11];

figure 
hold on
for i=1:numel(D)
    test=D{i};
    [ix,iy]=ind2sub(size(D),i);
    if isempty(test)
        list(ix,iy)={[]};
        spectrum_list(ix,iy)={[]};
    else
        
        factor=3.5;
        L=length(test(:,1));
        I=test(:,2);
        S_L=2^nextpow2(L);
        F_t=fft(I,S_L);
        dt=1/(abs(test(1,1)-test(2,1)));
        f=dt/2*linspace(0,1,S_L/2+1);
        p=F_t.*conj(F_t)/S_L;
        spectrum_list(ix,iy)={[f/f(end);(p(1:S_L/2+1)/p(S_L/2+1))']};
        
        if ix==1
            subplot(1,2,1)
            hold on
            plot(spectrum_list{ix,iy}(1,:),spectrum_list{ix,iy}(2,:)*base_p(1),'r')
        elseif ix==2
            subplot(1,2,2)
            hold on
            plot(spectrum_list{ix,iy}(1,:),spectrum_list{ix,iy}(2,:)*base_p(2),'g')
        elseif ix==3
            subplot(1,2,2)
            hold on
            plot(spectrum_list{ix,iy}(1,:),spectrum_list{ix,iy}(2,:)*base_p(3),'m')
        elseif ix==4
            subplot(1,2,1)
            hold on
            plot(spectrum_list{ix,iy}(1,:),spectrum_list{ix,iy}(2,:)*base_p(4),'b')
        elseif ix==5
            subplot(1,2,2)
            hold on
            plot(spectrum_list{ix,iy}(1,:),spectrum_list{ix,iy}(2,:)*base_p(5),'c')
        end
        
        ind=round((S_L/2+1)/factor);
        f(1:ind)=[];
        f(end-ind:end)=[];
        f=log(f);
        
        f_t=log(p(ind+1:S_L/2-ind))';
        A=zeros(length(f_t),2);
        A(:,1)=f';
        A(:,2)=1;
        b=f_t';
        [x,~]=lsqr(A,b,1e-12,500);
        list(ix,iy)={x(1)};
    end
end
hold off

fs=6;
low_list1=[list{1,:}];
low_list1=sort(low_list1);
ind1=round(length(low_list1)/fs);
low_list1(1:ind1)=[];
low_list1(end-ind1:end)=[];

low_list2=[list{4,:}];
low_list2=sort(low_list2);
ind2=round(length(low_list2)/fs);
low_list2(1:ind2)=[];
low_list2(end-ind2:end)=[];

high_list1=[list{2,:}];
high_list1=sort(high_list1);
ind3=round(length(high_list1)/fs);
high_list1(1:ind3)=[];
high_list1(end-ind3:end)=[];

high_list2=[list{3,:}];
high_list2=sort(high_list2);
ind4=round(length(high_list2)/fs);
high_list2(1:ind4)=[];
high_list2(end-ind4:end)=[];

high_list3=[list{5,:}];
high_list3=sort(high_list3);
ind5=round(length(high_list3)/fs);
high_list3(1:ind5)=[];
high_list3(end-ind5:end)=[];

mean_list=[mean(-low_list1),mean(-low_list2),mean(-high_list1-0.35),mean(-high_list2),mean(-high_list3)];
std_list=[var(-low_list1),var(-low_list2),var(-high_list1),var(-high_list2),var(-high_list3)];

figure
hold on
errorbar(2,mean_list(1),std_list(1))
errorbar(2.05,mean_list(2),std_list(2))
errorbar(2.1,mean_list(3),std_list(3))
errorbar(2.15,mean_list(4),std_list(4))
errorbar(2.2,mean_list(5),std_list(5))
% scatter(abs(low_list1),ones(length(low_list1),1)*2)
% scatter(abs(low_list2),ones(length(low_list2),1)*2.05)
% scatter(abs(high_list1+0.55),ones(length(high_list1),1)*2.1)
% scatter(abs(high_list2),ones(length(high_list2),1)*2.15)
% scatter(abs(high_list3),ones(length(high_list3),1)*2.2)
hold off
        
%% Fitting 2.1/2.15/2.2 V  
list_cell={};
load('C:\D\Program Files\Project\Matlab Project\Degradation\Overall\Results\2.1Vfit.mat')
list_cell{1,1}=[Result_list(:,1),Result_list(:,3)];
list_cell{1,2}=[D{chose_x,chose_y}(:,1),k];
load('C:\D\Program Files\Project\Matlab Project\Degradation\Overall\Results\2.15Vfit.mat')
list_cell{2,1}=[Result_list(:,1),Result_list(:,3)];
list_cell{2,2}=[D{chose_x,chose_y}(:,1),k];
load('C:\D\Program Files\Project\Matlab Project\Degradation\Overall\Results\2.2Vfit.mat')
list_cell{3,1}=[Result_list(:,1),Result_list(:,3)];
list_cell{3,2}=[D{chose_x,chose_y}(:,1),k];

figure
hold on
plot(list_cell{1,1}(:,1),list_cell{1,1}(:,2))
plot(list_cell{1,2}(:,1),list_cell{1,2}(:,2))
plot(list_cell{2,1}(:,1),list_cell{2,1}(:,2))
plot(list_cell{2,2}(:,1),list_cell{2,2}(:,2))
plot(list_cell{3,1}(:,1),list_cell{3,1}(:,2))
plot(list_cell{3,2}(:,1),list_cell{3,2}(:,2))
hold off


%% FR0709
s=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709.xlsx');
Data=s(:,1);
[p1,data1,final1,curve1]=Cluster_fitting(Data,[],[]);
Data=s(:,3);
[p2,data2,final2,curve2]=Cluster_fitting(Data,[],[]);
Data=s(:,5);
[p3,data3,final3,curve3]=Cluster_fitting(Data,4.4,[]);
Data=s(:,7);
[p4,data4,final4,curve4]=Cluster_fitting(Data,[],final3(1));
Data=s(:,9);
[p5,data5,final5,curve5]=Cluster_fitting(Data,4.4,[]);

figure
hold on
plot(curve1(1,:),curve1(2,:),'r')
plot(data1,p1,'*')
plot(curve2(1,:),curve2(2,:),'g')
plot(data2,p2,'+')
plot(curve3(1,:),curve3(2,:),'b')
plot(data3,p3,'*')
plot(curve4(1,:),curve4(2,:))
plot(data4,p4,'+')
plot(curve5(1,:),curve5(2,:),'y')
plot(data5,p5,'*')
hold off


%% FR2729
s=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR2729.xlsx');
Data=s(:,1);
[p1,data1,final1,curve1]=Cluster_fitting(Data,[],[]);
Data=s(:,3);
[p2,data2,final2,curve2]=Cluster_fitting(Data,[],[]);
Data=s(:,5);
[p3,data3,final3,curve3]=Cluster_fitting(Data,[],[]);
figure
hold on
plot(curve1(1,:),curve1(2,:),'r')
plot(data1,p1,'*')
plot(curve2(1,:),curve2(2,:),'g')
plot(data2,p2,'+')
plot(curve3(1,:),curve3(2,:),'b')
plot(data3,p3,'*')
hold off

%% FR0709 Bi-Model
s=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709.xlsx');
Data=s(:,3);
[p1,data1,final1,curve1]=Cluster_fitting_bicluster(Data,1/3,[],[],[]);
Data=s(:,1);
[p2,data2,final2,curve2]=Cluster_fitting_bicluster(Data,1/3,[final1(1),final1(2)],[final1(3),final1(4)],[]);
Data=s(:,9);
[p3,data3,final3,curve3]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s(:,13);
[p4,data4,final4,curve4]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s(:,5);
[p5,data5,final5,curve5]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);

figure
hold on
plot(curve1(1,:),curve1(2,:),'r')
plot(data1,p1,'*')
plot(curve2(1,:),curve2(2,:),'g')
plot(data2,p2,'+')
plot(curve3(1,:),curve3(2,:),'b')
plot(data3,p3,'*')
plot(curve4(1,:),curve4(2,:))
plot(data4,p4,'+')
plot(curve5(1,:),curve5(2,:),'y')
plot(data5,p5,'*')
hold off


%% FR2729 Bi-Model
s=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR2729.xlsx');
Data=s(:,1);
[p1,data1,final1,curve1]=Cluster_fitting_bicluster(Data,1/3,[0.0705862385463146,3.27097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s(:,3);
[p2,data2,final2,curve2]=Cluster_fitting_bicluster(Data,1/3,[0.0846748471866453,6.13485294039499],[4.72806830442574,2.33737021383377],[]);
Data=s(:,5);
[p3,data3,final3,curve3]=Cluster_fitting_bicluster(Data,1/3,[0.0705862385463146,3.27097516048798],[5.17466525324745,5.36075366710877],[]);

figure
hold on
plot(curve1(1,:),curve1(2,:),'r')
plot(data1,p1,'*')
plot(curve2(1,:),curve2(2,:),'g')
plot(data2,p2,'+')
plot(curve3(1,:),curve3(2,:),'b')
plot(data3,p3,'*')
hold off

%% FR2729 ！！FR0709 Bi-Model
s1=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR2729.xlsx');
s2=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709.xlsx');
% s3=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0103.xlsx');
s4=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR37.xlsx');
Data=s2(:,5);
[p1,data1,final1,curve1]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s2(:,9);
[p4,data4,final4,curve4]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s1(:,1);
[p2,data2,final2,curve2]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s1(:,9);
[p3,data3,final3,curve3]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
% Data=s3(:,1);
% [p6,data6,final6,curve6]=Cluster_fitting_bicluster(Data,1/3,[0.0705862385463146,3.27097516048798],[5.17466525324745,5.36075366710877],[]);
% Data=s3(:,5);
% [p5,data5,final5,curve5]=Cluster_fitting_bicluster(Data,1/3,[0.0705862385463146,3.27097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s4(:,1)*log(3.5);
[p7,data7,final7,curve7]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s4(:,5)*log(10.5);
[p8,data8,final8,curve8]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);

figure
hold on
plot(curve1(1,:),curve1(2,:),'r')
plot(data1,p1,'*')
plot(curve2(1,:),curve2(2,:),'g')
plot(data2,p2,'+')
plot(curve3(1,:),curve3(2,:),'b')
plot(data3,p3,'*')
plot(curve4(1,:),curve4(2,:),'m')
plot(data4,p4,'+')
% plot(curve5(1,:),curve5(2,:),'r')
% plot(data5,p5,'*')
% plot(curve6(1,:),curve6(2,:),'g')
% plot(data6,p6,'+')
plot(curve7(1,:),curve7(2,:),'y')
plot(data7,p7,'*')
plot(curve8(1,:),curve8(2,:),'c')
plot(data8,p8,'+')
hold off

%% FR0709 Bi-Model Thermal Effect
s=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709_position.xlsx');
Data=s(:,17);
[p1,data1,final1,curve1]=Cluster_fitting_bicluster(Data,1/3,[0.612689089419631,9.85145459807679],[2.72188447050496,14.8378784350425],[]);
Data=s(:,19);
[p2,data2,final2,curve2]=Cluster_fitting_bicluster(Data,1/3,[0.612689089419631,9.85145459807679],[2.72188447050496,14.8378784350425],[]);
Data=s(:,9);
[p3,data3,final3,curve3]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s(:,13);
[p4,data4,final4,curve4]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s(:,5);
[p5,data5,final5,curve5]=Cluster_fitting_bicluster(Data,1/3,[0.0825862385463146,2.87097516048798],[5.17466525324745,5.36075366710877],[]);
Data=s(:,21);
[p6,data6,final6,curve6]=Cluster_fitting_bicluster(Data,1/3,[0.612689089419631,9.85145459807679],[2.72188447050496,14.8378784350425],[]);

figure
hold on
plot(curve1(1,:),curve1(2,:),'r')
plot(data1,p1,'*')
plot(curve2(1,:),curve2(2,:),'g')
plot(data2,p2,'+')
plot(curve6(1,:),curve6(2,:),'b')
plot(data6,p6,'*')
plot(curve4(1,:),curve4(2,:),'c')
plot(data4,p4,'+')
plot(curve3(1,:),curve3(2,:),'m')
plot(data3,p3,'*')
plot(curve5(1,:),curve5(2,:),'y')
plot(data5,p5,'*')
hold off

%% FR2729 ！！FR0709 Histrogram
s1=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709_position.xlsx');
s2=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR2729_position.xlsx');
data=[];
FR0709=[];
for i=2:2:10
    data=s1(:,i);
    data=data(~isnan(data));
    FR0709=[FR0709;data];
    data=[];
end
FR2729=[];
for i=2:2:6
    data=s2(:,i);
    data=data(~isnan(data));
    FR2729=[FR2729;data];
    data=[];
end
[h1,n1]=hist(FR0709,20);
h1=h1/sum(h1);
[h2,n2]=hist(FR2729,20);
h2=h2/sum(h2);
h=[h1',h2'];
n=[n1',n2'];
bar(n,h);

%% FR2729 Thermal Histrogram
s1=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709_position.xlsx');
data=[];
FR0709=[];
for i=4:2:10
    data=s1(:,i);
    data=data(~isnan(data));
    FR0709=[FR0709;data];
    data=[];
end
FR0709T=[];
for i=18:2:24
    data=s1(:,i);
    data=data(~isnan(data));
    FR0709T=[FR0709T;data];
    data=[];
end
FR0709T=1-FR0709T;
FR0709=1-FR0709;
[h1,n1]=hist(FR0709,16);
h1=h1/sum(h1);
[h2,n2]=hist(FR0709T,16);
h2=h2/sum(h2);
h=[h1',h2'];
n=[n1',n2'];
bar(n,h);


%% TDDB current statistics

load('C:\D\Program Files\Project\Matlab Project\Degradation\Overall\Results\TDDB2.1V_100.mat')
Data=load('C:\D\Program Files\Project\Matlab Project\Clustering data processing\Fourier\current.mat');
D=Data.res;
x=5;
s_current=[];
s_time=[];
for i=1:length(D(x,:))
    sample=D{x,i};
    if isempty(sample)
        continue
    end
    s_current=[s_current,sample(2,2)];
    s_time=[s_time,sample(end,1)];
end


current_mean=mean(s_current);
time_mean=mean(s_time)*12;
current_std=std(s_current);
time_std=std(s_time)*12;

s_list=[];
figure
hold on
for i=1:length(final_list(:,1))
    factor_c=current_mean/2;
    factor_t=normrnd(time_mean,time_std/2.2);
    s=cell2mat(final_list(i,2));
    s(:,1)=s(:,1)/s(end,1)*factor_t;
    s(:,2)=s(:,2)/s(end,2)*factor_t;
    s(:,3)=s(:,3)/s(1,3)*factor_c;
    if s(:,3)~=2e-4
        s(end,3)=2e-4;
    end
    s_list=[s_list,s(end,1)];
    u1=[1;s(:,1)];
    u2=[s(1,3)/1e8;s(:,3)];
    plot(u1,u2,'b')
end

for i=1:length(D(x,:))
    sample=D{x,i};
    if isempty(sample)
        continue
    end
    sample(end,2)=2e-4;
    plot(sample(:,1)*12,sample(:,2),'r')
end
hold off

%% Simulation Clustering
load('C:\D\Program Files\Project\Matlab Project\Degradation\Overall\Results\TDDB2.1V_100.mat')
ss=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709_position.xlsx');

s_time=ss(:,5);
s_time=sort(s_time);
s_time=s_time(~isnan(s_time));

k=sort(ss(:,9));
k=k(~isnan(k));
kk=k(end-1);

s=cell2mat(final_list(:,1));
data=sort(s(:,1));

% ind=round(length(data)/3);
% ff=data(end-ind+1:1:end);
% ff=ff-data(end-ind);
% ff=ff*4/5;
% data(end-ind+1:1:end)=ff+data(end-ind);

o=data-data(1);
o=o/o(end)*(s_time(end)-s_time(1))+s_time(1);
data=o;

n=1:length(data);
p=(n-0.3)/(length(data)+0.4);
u=log(-log(1-p));
% plot(log(data),u,'*')

% [p4,data4,final4,curve4]=Cluster_fitting_bicluster(data,1/3,[],[],[]);
[ss,data,final,curve]=Least_fitting_square_clustering(data,1/3,[],[]);

figure
hold on
plot(curve(1,:)/curve(1,end)*kk,curve(2,:))
plot(data/data(end)*kk,u,'*')

hold off


%% N fitting

s1=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\T63.xlsx');
r=s1(2:3,1:5);
pp=s1(5:6,1:5);
r(2,1:2)=r(2,1:2)*1.2;
u=[2.0,2.05,2.1,2.15,2.2];
A=zeros(5,2);
A(:,1)=log(u');
A(:,2)=1;
b1=log(r(1,:)');
b2=log(r(2,:)');
[x1,~]=lsqr(A,b1,1e-12,500);
[x2,~]=lsqr(A,b2,1e-12,500);
i=1:5:2000;
f1=exp((log(i)-x1(2))/x1(1));
f2=exp((log(i)-x2(2))/x2(1));

%
% r1=s1(:,8:2:10);
% r2=s1(:,13:2:15);
% u1=[2.1,2.2];
% A1=ones(2,2);
% A1(:,1)=log(u1');
% bx1=log(r1(1,:)');
% bx2=log(r1(2,:)');
% [xx1,~]=lsqr(A1,bx1,1e-12,500);
% [xx2,~]=lsqr(A1,bx2,1e-12,500);
% fx1=exp((log(i)-xx1(2))/xx1(1));
% fx2=exp((log(i)-xx2(2))/xx2(1));
% 
% 
% A2=ones(2,2);
% A2(:,1)=log(u1');
% by1=log(r2(1,:)');
% by2=log(r2(2,:)');
% [xy1,~]=lsqr(A2,by1,1e-12,500);
% [xy2,~]=lsqr(A2,by2,1e-12,500);
% fy1=exp((log(i)-xy1(2))/xy1(1));
% fy2=exp((log(i)-xy2(2))/xy2(1));



figure 
hold on
% plot(f1,i,'r')
% plot(u,r(1,:),'*')
% plot(f2,i,'b')
% plot(u,r(2,:),'*')


plot(u,pp(1,:))
plot(u,pp(1,:),'*')
plot(u,pp(2,:))
plot(u,pp(2,:),'*')

% plot(u1,r1(1,:),'*')
% plot(fx1,i,'m')
% plot(u1,r1(2,:),'*')
% plot(fx2,i,'g')
% 
% plot(u1,r2(1,:),'*')
% plot(fy1,i,'c')
% plot(u1,r2(2,:),'*')
% plot(fy2,i,'k')
hold off

%% Process viaration
s1=xlsread('C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\T63.xlsx');
r=s1(:,18:19);
u=[1,2,3,4];
figure
hold on
plot(u,r(:,1))
plot(u,r(:,1),'*')
plot(u,r(:,2))
plot(u,r(:,2),'*')
hold off


a1=[1.2,0.43,0.32,0.08];
a2=[12.4,5.6,4.13,2.87];
b=[1,2,3,4];
figure
hold on
plot(b,a1)
plot(b,a1,'*')
plot(b,a2)
plot(b,a2,'*')
hold off

%% Initial R distribution and Time distribution

Data=load('C:\D\Program Files\Project\Matlab Project\Clustering data processing\Fourier\current.mat');
D=Data.res;
list_g=cell(size(D));
list_b=cell(size(D));
spectrum_list=cell(size(D));
for i=1:numel(D)
    test=D{i};
    [ix,iy]=ind2sub(size(D),i);
    if isempty(test)||length(test)<5
        list_g(ix,iy)={[]};
        list_b(ix,iy)={[]};
        spectrum_list(ix,iy)={[]};
    else
        if ix==1
            test(:,2)=2./test(:,2);
        elseif ix==2
            test(:,2)=2.1./test(:,2);
        elseif ix==3
            test(:,2)=2.2./test(:,2);
        elseif ix==4
            test(:,2)=2.05./test(:,2);
        elseif ix==5
            test(:,2)=2.15./test(:,2);
        end
        ind=round(length(test)/8);
        test(1:ind,:)=[];
        test(end-ind:end,:)=[];
        A=ones(length(test),2);
        A(:,1)=log(test(:,2));
        A(:,2)=1;
        b=log(test(:,1));
        [x,~]=lsqr(A,b,1e-12,1000);
        list_g(ix,iy)={x(1)};
        list_b(ix,iy)={exp(-x(2)/x(1))};
    end
end
fit_r=[mean([list_b{1,:}]), mean([list_b{4,:}]), mean([list_b{2,:}]), mean([list_b{5,:}]), mean([list_b{3,:}])];
A_r=ones(5,2);
A_r(:,1)=log(fit_r)';
b_r=[2,2.05,2.1,2.15,2.2]';
[x_r,~]=lsqr(A_r,b_r,1e-12,1000);
lv1=exp((1-x_r(2))/x_r(1));
figure
hold on
errorbar(2,mean([list_b{1,:}]),std([list_b{1,:}]))
errorbar(2.1,mean([list_b{2,:}]),std([list_b{2,:}]))
errorbar(2.2,mean([list_b{3,:}]),std([list_b{3,:}]))
errorbar(2.05,mean([list_b{4,:}]),std([list_b{4,:}]))
errorbar(2.15,mean([list_b{5,:}]),std([list_b{5,:}]))
hold off



