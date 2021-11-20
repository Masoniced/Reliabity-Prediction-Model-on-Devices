%% Resistor Network Solver
function [V,I,E,T_v,dT]=RECAL_i(VI,C_diag,s)
%tic
global P n1 n2 n3 aphla_C Location_lable Kth maxit Thermal_contribution Initial_T T_hk Map 

%% update of conductivity

%%
% f=0.2;
s_const=2e-4;
G=P'*C_diag*P;
% G(2,:)=[];
% G(:,2)=[];
L=ichol(G,struct('type','ict','droptol',1e-3,'diagcomp',500));
in=zeros(length(G),1);
lb=1;
in(1)=lb;
in(2)=-lb;
tol=max(min(VI/(n1*n2*n3*1000),min(s)/(n1*n2*n3)),1e-10);
r=1;
while r>=0.05
    [V,~,r]=minres(G,in,tol,maxit,L,L');
    if r>=0.05
        maxit=round(maxit+maxit/2);
    end
end
V=V-V(2);
factor=VI/V(1);
V=V*factor;
IB=-C_diag*P*V;
I=lb*factor;
IB(abs(IB)>I)=I;
V([1,2])=[];
V=reshape(V,[n1,n2,n3-2]);
V=V.*Map;
[gx,gy,gz]=gradient(V);
E=(abs(gx)+abs(gy)+abs(gz))./(3*aphla_C);
%% E smoothing

%E(:,:,1:T_hk)=E(:,:,1:T_hk).*Map;

%%
thermal=Thermal_contribution(1:end-1);
T_v=Initial_T+IB.^2./s.*Kth.*s_const.*thermal;
T_v=[T_v;0];
dT=zeros(n1,n2,n3-2);
for i=1:n1*n2*(n3-2)
    dT(i)=T_v(Location_lable(i,1))+T_v(Location_lable(i,2))+T_v(Location_lable(i,3))+T_v(Location_lable(i,4))+T_v(Location_lable(i,5))+T_v(Location_lable(i,6));
end
dT(1,1,:)=dT(1,1,:)/4;
dT(1,n2,:)=dT(1,n2,:)/4;
dT(n1,n2,:)=dT(n1,n2,:)/4;
dT(n1,1,:)=dT(n1,1,:)/4;
dT(1,2:n2-1,:)=dT(1,2:n2-1,:)/5;
dT(n1,2:n2-1,:)=dT(n1,2:n2-1,:)/5;
dT(2:n1-1,1,:)=dT(2:n1-1,1,:)/5;
dT(2:n1-1,n2,:)=dT(2:n1-1,n2,:)/5;
dT(2:n1-1,2:n2-1,:)=dT(2:n1-1,2:n2-1,:)/6;
%toc
end