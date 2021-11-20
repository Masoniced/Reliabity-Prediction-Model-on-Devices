%% MMC Degradation Module
global n1 n2 n3 VI
global Lable_mat label Label_record M
global P IVR_lable R_lable Location_lable location_conduct tunneling_lable Thermal_contribution R
global Ks p_h p_s kh ks_v K D_K aphla_C alph kh_gb ks_p
global Conduct E_con e
global Ea_h Ea_s Er_h Er_s Eh_h Eh_s 
global Kth Kth_HK Kth_IL
global ovcp_IL ovcp_HK Entro_HK Entro_IL
global maxit
global Initial_T T_hk P_matrix

tic
%% Specific parameter setting
% Defct distance
aphla_C=2e-10; % m
% Typical vibration frequency
Ks=1e13; % Hz
% Electron charge
e=1.6e-19; % C 
% Free activation of IL
Ea_s=0.7*e; % J
% Free activation of HK
Ea_h=1.2*e; % J
% DM of IL
p_s=4.275*e*1e-10; % e*`A
% DM of HK
p_h=10.2*e*1e-10; % e*`A
% Permitt of HK
kh=25;
% Permitt of HK GB
kh_gb=27;
% Permitt of virgin IL
ks_v=6;
% Permitt of Post BD IL
ks_p=8.5;
% Boltzmann constant
K=1.38e-23;
% Conduction
Conduct=[0.05*aphla_C,1.67*aphla_C]; % D. Berco1, T.Y. Tseng, J. Comput. Electron.,10.1007/s10825-015-0736-7;
% Conductive energy
E_con=[0.05*e,0.12*e]; % D. Berco1, T.Y. Tseng, J. Comput. Electron.,10.1007/s10825-015-0736-7;
%
alph=3.3e-10; % m
% Thermal resistance of HK
Kth_HK=1.5625e10/8e-10*aphla_C; % K/W (thermal conductivity)*(unit cell length)
% Thermal resistance of IL
Kth_IL=9.615e10/8e-10*aphla_C; % K/W (thermal conductivity)*(unit cell length)
%% 
% Er_h=0.95*e; % HK Recombination E
% Er_s=1.6*e; % IL Recombination E
% Eh_h=1.3*e; % HK Hopping E
% Eh_s=1.2*e; % IL Hopping E
Er_h=1.0*e; % HK Recombination E
Er_s=1.95*e; % IL Recombination E
Eh_h=1.1*e; % HK Hopping E
Eh_s=1.9*e; % IL Hopping E
%%
% Potential of HK OV
ovcp_HK=6.38*e; %J      %%  W. L. Scopel, Anto?nio J. R. da Silva, W. Orellana, and A. Fazzio,APL,84,9,2004
% Potential of IL OV
ovcp_IL=5.16*e; %J      %%  W. L. Scopel, Anto?nio J. R. da Silva, W. Orellana, and A. Fazzio,APL,84,9,2004
% Entropy of HK
Entro_HK=2.73*10^6*aphla_C^3; % Discussion of Theoretical Studies: Sections I-VI  Harold L. Schick 59.33 J/K*mol = 2.73 J/K*cm^3
% Entropy of IL
Entro_IL=1.59*10^6*aphla_C^3; 
%% Geometry Dimension Setting and Clustering Label Matrix (Changed)
n1=190;
n2=100;
n3=11;
T_hk=7;
M=zeros(n1,n2,n3-2); % Occupation matrix
Label_record={}; 
Lable_mat=zeros(size(M));
label=0;
P_matrix=ones(n1,n2,n3-2);

%% Initialise Basic Parameters for Device Environment
Initial_T=400;

C_hk=Conduct(1)*exp(-E_con(1)/(K*Initial_T));
C_il=Conduct(2)*exp(-E_con(2)/(K*Initial_T));
T=ones(n1,n2,n3-2)*Initial_T;
T_v=ones((3*n1*n2*n3-n1*n2-n2*n3-n1*n3-(2*n1*n2-n1-n2)*2)+1,1)*Initial_T;
T_v(end)=0;
tunneling_lable=zeros(2*n1-1,2*n2-1,2*n3-1);
OV_list=cell(1,2);
PIT_p=1/n2;

VI=2;

maxit=300;
%% 
% Initial network setting 
[P,C_diag,IVR_lable,R_lable,Location_lable,location_conduct,s,Kth,Thermal_contribution,D_K,Pre_check_value,R]=Initial_setting(n1,n2,n3,T_hk,C_hk,C_il,Kth_HK,Kth_IL,PIT_p);
toc
% Initial E/V/I calculation
[V,I,E,T_v,dT]=RECAL_i(VI,C_diag,s);
% Initial Probability setting
[Probability,total_rate]=Update(E,T,OV_list,T_hk);

%% Geometrical List
x=load('C:\Users\Mason\Desktop\P10.txt'); % Line from IL to HK
list=zeros(length(x(:,1)),2);
[list(:,1),ind_list]=sort(x(:,1));
tranist=x(:,4);
list(:,2)=tranist(ind_list);
[~,ind_list]=max(abs(diff(diff(list(:,2))./diff(list(:,1)))));
IL=list(1:ind_list+1,:);
HK=list(ind_list+2:end,:);
IL(:,2)=(IL(:,2)-IL(1,2))/(IL(end,2)-IL(1,2));
IL(:,1)=(IL(:,1)-IL(1,1))/(IL(end,1)-IL(1,1));
HK(:,2)=(HK(:,2)-HK(1,2))/(HK(end,2)-HK(1,2));
HK(:,1)=(HK(:,1)-HK(1,1))/(HK(end,1)-HK(1,1));
list_p=zeros(n3-2,1);
for i1=1:n3-2
    if i1<=T_hk
        asump=(T_hk+1-i1)/(T_hk+1);
        [~,ind_hk]=min(abs(HK(:,1)-asump));
        if HK(ind_hk,1)>=asump
            x2=HK(ind_hk,1);
            y2=HK(ind_hk,2);
            x1=HK(ind_hk-1,1);
            y1=HK(ind_hk-1,2);
        else
            x1=HK(ind_hk,1);
            y1=HK(ind_hk,2);
            x2=HK(ind_hk+1,1);
            y2=HK(ind_hk+1,2);
        end
        list_p(i1)=((asump-x1)/(x2-x1)*(y2-y1)+y1)/asump;
    else
        if i1==T_hk+1
            list_p(i1)=1;
        else
            asump=(n3-2-i1+1)/(n3-2-T_hk+1);
            [~,ind_il]=min(abs(IL(:,1)-asump));
            if IL(ind_hk,1)>=asump
                x2=IL(ind_il,1);
                y2=IL(ind_il,2);
                x1=IL(ind_il-1,1);
                y1=IL(ind_il-1,2);
            else
                x1=IL(ind_il,1);
                y1=IL(ind_il,2);
                x2=IL(ind_il+1,1);
                y2=IL(ind_il+1,2);
            end
            list_p(i1)=((asump-x1)/(x2-x1)*(y2-y1)+y1)/asump;  
        end
    end
    P_matrix(:,:,i1)=list_p(i1);
end