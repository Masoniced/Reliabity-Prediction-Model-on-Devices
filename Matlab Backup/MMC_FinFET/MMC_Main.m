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
global Initial_T

tic
%% Specific parameter setting
% Defct distance
aphla_C=4e-10; % m
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
n1=175;
n2=60;
n3=10;
T_hk=4;
M=zeros(n1,n2,n3-2); % Occupation matrix
Label_record={}; 
Lable_mat=zeros(size(M));
label=0;
%% Initialise Basic Parameters for Device Environment
Initial_T=400;

C_hk=Conduct(1)*exp(-E_con(1)/(K*Initial_T));
C_il=Conduct(2)*exp(-E_con(2)/(K*Initial_T));
T=ones(n1,n2,n3-2)*Initial_T;
T_v=ones((3*n1*n2*n3-n1*n2-n2*n3-n1*n3-(2*n1*n2-n1-n2)*2)+1,1)*Initial_T;
T_v(end)=0;
tunneling_lable=zeros(2*n1-1,2*n2-1,2*n3-1);
OV_list=cell(1,2);
PIT_p=1.5/n2;

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

%% MMC Experiments

criteria=1;
counter=0;
time=0;
mini_batch=10;
total_time=0;
converge_count=0;
I_list=[];
Emax_list=[];
I_list=[I_list,I];
[E_value,E_position]=max(E(:));
Emax_list=[Emax_list;E_value,E_position];
Check_Volum=[];

while criteria==1
    rng('shuffle')
    if counter==n1*n2*(n3-2)
        break
    end    
    seed=rand(mini_batch,3);
    delt_time=zeros(mini_batch,1);
    defect_configure=cell(mini_batch,3);
    batch_time=0;
    %% Update of Probability Vector
    for i=1:mini_batch
        fix=1;
        while fix==1
            if seed(i,1)<=total_rate(1,5)
                pick_y=1;
                defect_configure{i,3}='H';
            else
                pick_y=2;
                defect_configure{i,3}='L';
            end

            [~,pick_x]=min(abs(total_rate(:,pick_y+2)-seed(i,2)));
            process_Probability=Probability{pick_x,pick_y};
            if ~isempty(process_Probability)
                fix=0;
            else
                seed(i,:)=rand(1,3);
            end
        end
        [~,ind]=min(abs(process_Probability(:,1)-seed(i,3)));
        
        defect_configure{i,1}=pick_x;
        
        if pick_x==1
            [new_Probability]=Probability_vector_process(process_Probability,ind);

            if pick_y==1
                defect_configure{i,2}=Probability{pick_x,pick_y}(ind,2);
            else
                defect_configure{i,2}=Probability{pick_x,pick_y}(ind,2)+n1*n2*T_hk;
            end
            Probability{pick_x,pick_y}=new_Probability;
            
        elseif pick_x==2
            [new_Probability1]=Probability_vector_process(process_Probability,ind);
            
            process_Probability=Probability{pick_x+1,pick_y};
            [new_Probability2]=Probability_vector_process(process_Probability,ind);
            if pick_y==1
                defect_configure{i,2}=Probability{pick_x,pick_y}(ind,2);
            else
                defect_configure{i,2}=Probability{pick_x,pick_y}(ind,2)+n1*n2*T_hk;
            end
            
            Probability{pick_x,pick_y}=new_Probability1;
            Probability{pick_x+1,pick_y}=new_Probability2;
        else
            [new_Probability1]=Probability_vector_process(process_Probability,ind);
            
            process_Probability=Probability{pick_x-1,pick_y};
            [new_Probability2]=Probability_vector_process(process_Probability,ind);
            
            if pick_y==1
                defect_configure{i,2}(1)=Probability{pick_x,pick_y}(ind,2);
            else
                defect_configure{i,2}(1)=Probability{pick_x,pick_y}(ind,2)+n1*n2*T_hk;
            end
            Probability{pick_x,pick_y}=new_Probability1;            
            Probability{pick_x-1,pick_y}=new_Probability2;
            
            h_jp=zeros(1,3);
            h_jp(1)=rem(defect_configure{i,2}(1),n1);
            if h_jp(1)==0
                h_jp(1)=n1;
            end
            h_jp(2)=rem((defect_configure{i,2}(1)-h_jp(1))/n1,n2)+1;
            h_jp(3)=(defect_configure{i,2}(1)-h_jp(1)-n1*(h_jp(2)-1))/(n1*n2)+1;
            
            [destination]=Decide_hopping(h_jp);
            defect_configure{i,2}(2)=destination;
            if destination<=n1*n2*T_hk
                defect_configure{i,3}(2)='H';
            else
                defect_configure{i,3}(2)='L';
            end
                 
        end
        delt_time(i)=(log(1/seed(1))+log(1/seed(2))+log(1/seed(3)))*1/(total_rate(1,1)+total_rate(1,2)+total_rate(2,1)+total_rate(2,2)+total_rate(3,1)+total_rate(3,2));
        
        batch_time=batch_time+delt_time(i);
        
    end
    %% Decide Entropy for diverse events
    
    [V,temp_I,temp_E,T_v,T,OV_list,temp_Probability,temp_total_rate,criteria_index,fb_list,check_volum_f]=MMC_decide(defect_configure,mini_batch,VI,V,T_v,T,OV_list,T_hk);
    if ~isempty(temp_I)
        I=temp_I;
        E=temp_E;
        Probability=temp_Probability;
        total_rate=temp_total_rate;
        total_time=total_time+batch_time;
        if criteria_index==1
            criteria=0;
            Check_Volum=check_volum_f;
        end
        converge_count=0; % Set 0 for not converge
        counter=counter+1;
        I_list=[I_list,I];
        [E_value,E_position]=max(E(:));
        Emax_list=[Emax_list;E_value,E_position];
    else
        [Probability,total_rate]=Update(E,T,OV_list,T_hk);
        converge_count=converge_count+1;
        if converge_count>=5; % continous 5 minibatch as converge criteria
            criteria=0;
            fprintf('\n Convergence Approaching')
        end
        
    end
end
    
    
    
    
    
    