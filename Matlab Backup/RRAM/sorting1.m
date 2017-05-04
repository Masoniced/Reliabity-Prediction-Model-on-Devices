function [probability,podes,sped,Incretemp,T_rate,J_rate]=sorting1(position,P_time,C_time)
global clasc C O n1 T Cur Unitcon Res MD_ns MD_s Cp_ns Cp_s tk_ns tk_s tcr Initial_T 
%% Parameters setting
p=position; 
S_rate=zeros(1,17);
podes=zeros(17,3);
sped=zeros(1,17);
R_rate=zeros(17,4);
C_temp=T(p(1),p(2),p(3));
%% Concentration at intial position
if p(1)<=2||p(1)>=4*n1||p(2)<=2||p(2)>=4*n1
    probability=0;podes=[0,1,1];sped=0;Incretemp=0;
    return
else
    if (p(3)==2)&&(~(p(1)<3||p(1)>4*n1-1||p(2)<3||p(2)>4*n1-1))
        sample1=O(p(1)-2:1:p(1)+2,p(2)-2:1:p(2)+2,p(3)-1:1:p(3)+3);
        sample2=C(p(1)-2:1:p(1)+2,p(2)-2:1:p(2)+2,p(3)-1:1:p(3)+3);
        q1=length(find(sample1==3));
        q2=length(find(sample2==2));
        poC=q1/18;ratio=q1/(q1+q2);
    elseif (p(3)==1)&&(~(p(1)<3||p(1)>4*n1-1||p(2)<3||p(2)>4*n1-1))
        sample1=O(p(1)-2:1:p(1)+2,p(2)-2:1:p(2)+2,p(3):1:p(3)+4);
        sample2=C(p(1)-2:1:p(1)+2,p(2)-2:1:p(2)+2,p(3):1:p(3)+4);
        q1=length(find(sample1==3));
        q2=length(find(sample2==2));
        poC=q1/18;ratio=q1/(q1+q2);
    else
        sample1=O(p(1)-2:1:p(1)+2,p(2)-2:1:p(2)+2,p(3)-2:1:p(3)+2);
        sample2=C(p(1)-2:1:p(1)+2,p(2)-2:1:p(2)+2,p(3)-2:1:p(3)+2);
        q1=length(find(sample1==3));
        q2=length(find(sample2==2));
        poC=q1/18;ratio=q1/(q1+q2);
    end
end
%% Thermal update 
Idensity=(Cur/(2*pi*Unitcon^2))*(1/(sqrt((2*(p(1)-2*n1)^2+2*(p(2)-2*n1)^2+(p(3))^2))+sqrt(1/2)*(0.2*n1+1))^2);% changed by position
% Second order derivative of T
% Boundary thermal condition
if (p(1)<3||p(1)>4*n1-1||p(2)<3||p(2)>4*n1-1)
    sod_t_time=0;
% Normal thermal condition
else
    sod_t_time=12*Idensity^2*Res*(1+tcr*(C_temp-Initial_T))*(1/((sqrt((2*(p(1)-2*n1)^2+2*(p(2)-2*n1)^2+(p(3))^2))+sqrt(1/2)*(0.2*n1+1))^2*Unitcon^2))/((MD_ns)*(Cp_ns));
end
Incretemp=(Idensity^2*Res*(1+tcr*(C_temp-Initial_T))*(P_time)+sod_t_time*tk_ns*(P_time)^(2.5)/2)/((MD_ns)*(Cp_ns));
if Incretemp<0
    fprintf('Temperature error');
    fprintf('%d\n',p,P_time,C_time,sod_t);
    return
end
% First order derivative of T
% Boundary thermal condition
th_vector=-(4*Idensity^2*Res*(1+tcr*(C_temp-Initial_T))*(P_time))/((MD_ns)*(Cp_ns)*(sqrt((2*(p(1)-2*n1)^2+2*(p(2)-2*n1)^2+(p(3))^2))+sqrt(1/2)*(0.2*n1+1)))-(6*sod_t_time*tk_ns*(P_time)^(2.5)/2)/((MD_ns)*(Cp_ns)*(sqrt((2*(p(1)-2*n1)^2+2*(p(2)-2*n1)^2+(p(3))^2))+sqrt(1/2)*(0.2*n1+1)));
fod_tx=th_vector*(1/(sqrt((2*(p(1)-2*n1)^2+2*(p(2)-2*n1)^2+(p(3))^2))+sqrt(1/2)*(0.2*n1+1)))*(sqrt(2)*(p(1)-2*n1));
fod_ty=th_vector*(1/(sqrt((2*(p(1)-2*n1)^2+2*(p(2)-2*n1)^2+(p(3))^2))+sqrt(1/2)*(0.2*n1+1)))*(sqrt(2)*(p(2)-2*n1));
fod_tz=th_vector*(1/(sqrt((2*(p(1)-2*n1)^2+2*(p(2)-2*n1)^2+(p(3))^2))+sqrt(1/2)*(0.2*n1+1)))*(p(3)+sqrt(1/2)*(0.2*n1+1));
th_direction=[fod_tx,fod_ty,fod_tz];
%% Decide the normalized total rates
for i=1:17;
    [podes(i,:),S_rate(i),sped(i),R_rate(i,:)]=searching5(p,clasc(i,:),poC,th_direction,Idensity,Incretemp);
end
a=find(podes==0);
if ~isempty(a)
    if size(a)==17;
        podes=zeros(1,3);
        J_rate=zeros(1,4);
        sped=0;
        probability=0;
        fprintf('Temperature error');
        fprintf('%d\n',p,podes,C_time,sod_t);
    else
        podes(a,:)=[];
        R_rate(a,:)=[];
        S_rate(a)=[];
        sped(a)=[];
        J_rate=[sum(R_rate(:,1)),sum(R_rate(:,2)),sum(R_rate(:,3)),sum(R_rate(:,4))];
        T_rate=sum(S_rate);
        N_rate=S_rate/T_rate;
        t=length(N_rate);
        probability=zeros(1,t);
        probability(1)=N_rate(1);
        for i=2:t
            probability(i)=probability(i-1)+N_rate(i);
        end
    end
else
    J_rate=[sum(R_rate(:,1)),sum(R_rate(:,2)),sum(R_rate(:,3)),sum(R_rate(:,4))];
    T_rate=sum(S_rate);
    N_rate=S_rate/T_rate;
    t=length(N_rate);
    probability=zeros(1,t);
    probability(1)=N_rate(1);
    for i=2:t
        probability(i)=probability(i-1)+N_rate(i);
    end
end
end
