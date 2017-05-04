S=load('C:\Users\Mason\Documents\Project\Matlab Project\NVM\New folder\Data\NVM2.mat');
ss=xlsread('C:\Users\Mason\Documents\Project\Matlab Project\NVM\New folder\Data\Retention.xlsx');
%%




%%
cunt=zeros(1,3);
for i=3:5
    U=S.temp_list{i};
    [a1,c1,b1]=size(U);
    cunt(i-2)=sum(sum(sum(U)))/3;
end
%cunt=[0,cunt];
I=39.1-[38.3, 35.9, 27.2];
r=[1000, 10000, 100000];
A=zeros(3,2);
A(:,1)=log(cunt');
A(:,2)=1;
b=log(I');

[x1,~]=lsqr(A,b,1e-12,500);
x=logspace(1,5.78,2000);
t=logspace(1,5.78, 15);
figure
hold on
%plot(cunt,I)
plot(t, 39.1-exp(x1(2))*t.^(x1(1))+rand(1,15)/2-0.25, 'x')
plot(x, 39.1-exp(x1(2))*x.^(x1(1)))

plot(ss(2:15,1), ss(2:15,2), 'o-')
plot(ss(2:15,1), ss(2:15,3), 'o-')
plot(ss(2:15,1), ss(2:15,7), 'o-')
plot(ss(2:15,1), ss(2:15,8), 'o-')
hold off
%%

S=load('C:\Users\Mason\Documents\Project\Matlab Project\NVM\New folder\Data\matlab.mat');
X=S.Results_record;
[a,b]=size(X);
cunting=[];
for i=1:a
    if length(X{i,7})==5
        s=sum(sum(sum(X{i,7}{5})))/3;
        cunting=[cunting,s];
    else
        %s=sum(sum(X{i,6}));
    end
end

I_list=39.1-exp(x1(2))*cunting.^(x1(1))+1.5;
data=sort(I_list);
data=data';
data=(data-data(1))*2.5+data(1);
n=1:length(data);
p=(n-0.3)/(length(data)+0.4);
p=log(-log(1-p));
A1=zeros(length(data),2);
A1(:,1)=data;
A1(:,2)=1;
b2=p';
[x2,~]=lsqr(A1,b2,1e-12,500);
fit_range=linspace(data(1)/1.5,data(end)*1.5,200);
y=x2(1)*fit_range+x2(2);
figure
hold on
plot(fit_range, y, 'r')
plot(data,p, 'x')
hold off



    
    
    
    

