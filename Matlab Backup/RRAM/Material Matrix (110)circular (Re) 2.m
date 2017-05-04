syms samp
global C I Diffusion O n1 n2 K T e Z Cur Qs clasc Da omega Module_S N Unitcon Res MD_ns MD_s Cp_ns Cp_s tk_ns tk_s tcr Initial_T pr th_ex Ks
%% Build up the structure along (110) direction and initial diffusion matrix (n1/10)*(n1/10)
c1(:,:,1)=[2,0,2,0,2;0,0,0,0,0;2,0,2,0,2;0,0,0,0,0;2,0,2,0,2];
c1(:,:,2)=[0,0,0,0,0;2,0,2,0,2;0,0,0,0,0;2,0,2,0,2;0,0,0,0,0];
c1(:,:,3)=[0,0,0,0,0;0,2,0,2,0;0,0,0,0,0;0,2,0,2,0;0,0,0,0,0];
c1(:,:,4)=[0,2,0,2,0;0,0,0,0,0;0,2,0,2,0;0,0,0,0,0;0,2,0,2,0];
n1=100;n2=201;%The size of structure is (4*n1+1)*(4*n1+1)*(2*n2-1)
C=repmat(c1,[n1 n1 n1+1]);
n3=5:5:5*n1-5;
C(n3,:,:)=[];
C(:,n3,:)=[];
C(:,:,4*n1+4)=[];
C(:,:,4*n1+3)=[];
C(:,:,4*n1+2)=[];%Si srtucture matrix
i1(:,:,1)=[1,3,1,3,1;3,1,3,1,3;1,3,1,3,1;3,1,3,1,3;1,3,1,3,1];
i1(:,:,2)=[3,1,3,1,3;1,3,1,3,1;3,1,3,1,3;1,3,1,3,1;3,1,3,1,3];
i1(:,:,3)=[1,3,1,3,1;3,1,3,1,3;1,3,1,3,1;3,1,3,1,3;1,3,1,3,1];
for i=4:(2*n2-1)
    i1(:,:,i)=[1,1,1,1,1;1,1,1,1,1;1,1,1,1,1;1,1,1,1,1;1,1,1,1,1];
end
[xx,yy]=meshgrid(1:4*n1+1);
Si=((xx-2*n1-1).^2+(yy-2*n1-1).^2)<=(0.2*n1)^2/2;
Si=Si-ones(4*n1+1,4*n1+1);
Sind=sub2ind(size(Si),find(Si==-1));
I=repmat(i1,[n1 n1]);
I(n3,:,:)=[];I(:,n3,:)=[];
s1=I(:,:,1);s1(Sind)=1;I(:,:,1)=s1;
s2=I(:,:,2);s2(Sind)=1;I(:,:,2)=s2;
s3=I(:,:,3);s3(Sind)=1;I(:,:,3)=s3;%Intitial srtucture matrix
O=I;%Ocuppation matrix
T=zeros(4*n1+1,4*n1+1,2*n2-1)+300;
%% Basic parameters setting
% Unit constant length
Unitcon=(5.43e-10)/4; % m
% Boltzmann constant
K=1.38e-23;
% Initial temperature
Initial_T=300; % K
% Electron charge
e=1.6e-19; % C
% Effective charge for Ni
Z=0.62;
% Compliance current
Cur=5e-6;% A/m^2
% Intrinsic diffusion coefficient
Da=2e-7;% m^2/s
% Atomic volume
omega=1.21e-5; % m^3
% Young's module of NiSi
Module_S=1.6e11; % Pa=N/m
% Poisson ratio of NiSi
pr=0.42;
% Thermal expansion coefficient of NiSi
th_ex=1.6e-5; % /K
% Avogadro constant
N=6.02e23;
% Resistivity of NiSi
Res=3.4e-7; % m/S
% Temperature coefficient of resistivity of NiSi
tcr=0.001; % /K
% Density of NiSi
MD_ns=7200; %Kg/m^3
% Density of Si
MD_s=2329; %Kg/m^3
% Heat compacity of NiSi
Cp_ns=191.739;
% Heat compacity of Si
Cp_s=700;
% Heat conductivity of NiSi
tk_ns=33;
% Heat conductivity of Si
tk_s=149;
% Heat of transfer
Qs=0.2*e;
% Typical lattice frequency of NiSi
Ks=2.95e13;
%Diffusion 
Diffusion(1).dir=[1,1,0;1,-1,0;-1,1,0;-1,-1,0;0,0,2];Diffusion(1).spa=1/8;Diffusion(1).dis=2;Diffusion(1).aEnergy=0.76;%(100)
Diffusion(2).dir=[2,0,0;-2,0,0;0,2,0;0,-2,0;1,1,2;-1,1,2;1,-1,2;-1,-1,2];Diffusion(2).spa=sqrt(2)/4;Diffusion(2).dis=2*sqrt(2);Diffusion(2).aEnergy=0.26;%(110)
Diffusion(3).dir=[2,0,2;-2,0,2;0,2,2;0,-2,2];Diffusion(3).spa=sqrt(3)/4;Diffusion(3).dis=2*sqrt(3);Diffusion(3).aEnergy=0.425;%(111) 
Diffusion(4).sur=[];Diffusion(4).aEnergy=0.26;%Dependence of direction
clasc(1,:)=[1,1];clasc(2,:)=[1,2];clasc(3,:)=[1,3];clasc(4,:)=[1,4];clasc(5,:)=[1,5];% Sequence of sampling direction
clasc(6,:)=[2,1];clasc(7,:)=[2,2];clasc(8,:)=[2,3];clasc(9,:)=[2,4];clasc(10,:)=[2,5];clasc(11,:)=[2,6];clasc(12,:)=[2,7];clasc(13,:)=[2,8];
clasc(14,:)=[3,1];clasc(15,:)=[3,2];clasc(16,:)=[3,3];clasc(17,:)=[3,4];
%% Kinetic Monte Carlo
criteria=0;
C_time=0;
C_time1=0;
pr_in=1;
P_time=1e-18;
number1=0;
pos_t=[];
Para_record=[];
Shape_record.plane_view=[];
Shape_record.top_view=[];
Shape_record.angle=[];
while criteria~=1
    [x,y,z]=ind2sub(size(O),find(O==3));% Index of all initial positions
    rng('shuffle')
    k1=length(x);
    rng('shuffle')
    r1=rand(1,k1);
    r2=rand(1,k1);
    [p1,p2]=sort(r2);
    decide_fill=zeros(k1,3);
    t=zeros(1,k1);
    t_time=0;
    S=zeros(4*n1+1,4*n1+1,2*n2-1);
    S_rate=zeros(4*n1+1,4*n1+1,2*n2-1)+inf;
    add_temp=zeros(4*n1+1,4*n1+1,2*n2-1);
    K_rate=zeros(1,4);
    for ok=1:k1
        l=p2(ok);
       [probability,posdes,sped,Incretemp,T_rate,J_rate]=sorting2([x(l),y(l),z(l)],P_time,C_time);
       if (posdes(1)==0)
           break;
       end
       [min_value,index]=min(abs(probability-r1(l)));
       o=length(index);
       if o==1
           if probability(index)>r1(l)
               con_index=index;
           else
               con_index=index+1;
           end
       else
           if probability(index(o))<r1(l)
               con_index=index(o)+1;
           else
               dr=rand(1,1);subindex=1/o:1/o:1;
               [pass_para,pionts]=min(abs(subindex-dr));
               if subindex(points)>dr
                   con_index=index(points);
               else
                   con_index=index(points)+1;
               end
           end
       end
       decide_fill(ok,:)=posdes(con_index,:);
       S(decide_fill(ok,1),decide_fill(ok,2),decide_fill(ok,3))=2;
       t(ok)=(1/((T_rate*(2*omega/(N))^(2/3))*(3/(4*pi))^(2/3)*4*pi*2*(2*n2+1)^3))*log(1/r1(l));
       if S_rate(decide_fill(ok,1),decide_fill(ok,2),decide_fill(ok,3))==inf
           S_rate(decide_fill(ok,1),decide_fill(ok,2),decide_fill(ok,3))=1/t(ok);
       else
           S_rate(decide_fill(ok,1),decide_fill(ok,2),decide_fill(ok,3))=S_rate(decide_fill(ok,1),decide_fill(ok,2),decide_fill(ok,3))+1/t(ok);
       end
       t_time=t_time+t(ok);
       P_time=(t_time/ok)*k1*(pr_in)^2;% even distribution random parameter estimation
       K_rate=K_rate+J_rate;
       add_temp(x(l),y(l),z(l))=Incretemp;
    end
    add_temp=smothing(add_temp,x,y,z);
    el=length(find(S));
    pr_in=el/k1;
    S_rate=1./S_rate;
    number1=number1+1;
    t_time1=sum(sum(sum(S_rate)));
    totalr=sum(K_rate);
    % pos_t(number1)=t_time;
    pos_t1(number1)=t_time1;
    % C_time=C_time+t_time;
    C_time1=C_time1+t_time1;
    C_time=C_time1;
    O=O+S;
    T=T+add_temp;
    max_temp=max(max(max(T)));
    [Shape_record.angle(number1,1),Shape_record.angle(number1,2),Shape_record.plane_view(:,:,number1),Shape_record.top_view(:,:,number1)]=painting_record(O);
    Para_record(number1,1)=C_time;
    Para_record(number1,2)=K_rate(1)/totalr;
    Para_record(number1,3)=K_rate(2)/totalr;
    Para_record(number1,4)=K_rate(3)/totalr;
    Para_record(number1,5)=K_rate(4)/totalr;
    Para_record(number1,6)=totalr;
    Para_record(number1,7)=max_temp;
    fprintf('Total Time,percentage');
    fprintf('%d\n',C_time,C_time/5e-12,max_temp);
    clear r1 x y z t
    if ~isempty(find(O(:,:,2*n2-3)==3, 1))
        criteria=1;
    elseif C_time>5e-10
        criteria=1;
    elseif max_temp>450
        criteria=1;
    else
        criteria=0;
    end
end

