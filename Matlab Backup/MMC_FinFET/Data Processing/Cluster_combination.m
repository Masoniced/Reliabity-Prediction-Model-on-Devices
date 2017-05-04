function [t_G1,t_G2,t_G3,t_G4,t_G5,t_G6,F,x0,data,p_s]=Cluster_combination(Data,ratio,batch_n)
syms a1 b1 x a2 b2 y t
CDF=1-((1+1/a1*(t/x)^b1)^(-a1))*((1+1/a2*(t/y)^b2)^(-a2));
sym1=strcat('tor',num2str(2*batch_n-1));
sym2=strcat('tor',num2str(2*batch_n));
syms(strcat('tor',num2str(2*batch_n-1)));
syms(strcat('tor',num2str(2*batch_n)));
CDF=subs(CDF,x,sym1);
CDF=subs(CDF,y,sym2);
% s=xlsread('C:\D\Program Files\MATLAB\projeects\Clustering data processing\FinFET\125 TDDB FR09.xlsx');
% Data=s(:,2);
data=Data;
data=data(~isnan(data));
data=sort(data);
x0=[data(round(length(data)*0.63*ratio)),data(round(length(data)*0.63*(1-ratio)+length(data)*ratio))];

n=1:length(data);
p=(n-0.3)/(length(data)+0.4);
p_s=log(-log(1-p));     

fa1=diff(CDF,a1);
fb1=diff(CDF,b1);
ftor1=diff(CDF,sym1);
fa2=diff(CDF,a2);
fb2=diff(CDF,b2);
ftor2=diff(CDF,sym2);
Lf=matlabFunction(CDF);
G1=matlabFunction(fa1);
G2=matlabFunction(fa2);
G3=matlabFunction(fb1);
G4=matlabFunction(fb2);
G5=matlabFunction(ftor1);
G6=matlabFunction(ftor2);
t_G1=0;
t_G2=0;
t_G3=0;
t_G4=0;
t_G5=0;
t_G6=0;
F=0;


for i=1:length(data)
    t_G1=2*subs(G1,t,data(i)).*(subs(Lf,t,data(i))-p(i))+t_G1;
    t_G2=2*subs(G2,t,data(i)).*(subs(Lf,t,data(i))-p(i))+t_G2;
    t_G3=2*subs(G3,t,data(i)).*(subs(Lf,t,data(i))-p(i))+t_G3;
    t_G4=2*subs(G4,t,data(i)).*(subs(Lf,t,data(i))-p(i))+t_G4;
    t_G5=2*subs(G5,t,data(i)).*(subs(Lf,t,data(i))-p(i))+t_G5;
    t_G6=2*subs(G6,t,data(i)).*(subs(Lf,t,data(i))-p(i))+t_G6;     
    F=(subs(Lf,t,data(i))-p(i)).^2+F;
end

end

