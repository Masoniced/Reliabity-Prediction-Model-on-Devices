%% Test
clear
clc
load('C:\Users\Mason\Documents\Project\Matlab Project\NC\Pyproduct\initial3_70_70_28.mat')
% Data=load('C:\D\Program Files\Project\Matlab Project\Clustering data processing\Fourier\current.mat');
% D=Data.res;
% chose_x=5;
% chose_y=3;
% scale=D{chose_x,chose_y}(end-1,2)/D{chose_x,chose_y}(1,2);
% I_factor=D{chose_x,chose_y}(end-2,2);
% t_factor=D{chose_x,chose_y}(end,1);

%%
%% V distribution
ss=load('C:\Users\Mason\Desktop\Plot_2_0_data2.csv');
ss(:,2)=[];
ss(ss(:,1)>14,:)=[];
ss(ss(:,1)<-14,:)=[];
ss(ss(:,2)>7,:)=[];
ss(ss(:,2)<-3.74,:)=[];
ss(:,2)=ss(:,2)-2;
Map_list=zeros(n1,n3-2);
Map=zeros(n1,n2,n3-2);

for i=1:n1*(n3-2)
    
    [x,y]=ind2sub([n1,n3-2],i);
    if y<=T_hk
        u=((x-0.5)*28/n1-14-ss(:,1)).^2+((T_hk-y-0.5)*10/(T_hk)-5-ss(:,2)).^2;
    else
        u=((x-0.5)*28/n1-14-ss(:,1)).^2+(-(y-T_hk-0.5)*9.6/(T_hk)-5.4-ss(:,2)).^2;
    end
    [v,ind]=min(u);
    Map_list(x,y)=ss(ind,3);
    
end
U1=Map_list(:,1:T_hk);
U2=Map_list(:,T_hk+1:end);
HK_max=max(U1(:));
HK_min=min(U1(:));
IL_max=max(U2(:));
IL_min=min(U2(:));

for i=1:n1*(n3-2)
    [x,y]=ind2sub([n1,n3-2],i);
    if y<=T_hk
        Map_list(x,y)=Map_list(x,y)/((HK_max-(y/T_hk)*(HK_max-HK_min))*2);
    else
        Map_list(x,y)=Map_list(x,y)/((IL_max-((y-T_hk)/(n3-2-T_hk))*(IL_max-IL_min))*2);
    end
    Map(x,:,y)=Map_list(x,y);
end
%%

VI=4;

n_w=2.8e-8;
n_l=2.8e-7;
n_f=1;
t_o=1e-8;

v_factor_1=-56.0564;
v_factor_2=45.5015;


x_r=[32.8248319292306;-37.0387213298481];
slope=0.8588;
Ipre_mean_variation_factor=0.2448;

Initial_Rp=(VI/(VI^(x_r(1))*exp(x_r(2))))*(n_w*n_l*n_f)/t_o;
Initial_Ipre_mean=(VI^(x_r(1))*exp(x_r(2)));
viaration_Ipre=Initial_Ipre_mean*Ipre_mean_variation_factor;
Initial_Ipre=normrnd(Initial_Ipre_mean,viaration_Ipre);


Ea_s=0.7*e; % J
Er_h=1.0*e; % HK Recombination E
Er_s=1.95*e; % IL Recombination E
Eh_h=1.1*e; % HK Hopping E
Eh_s=1.9*e; % IL Hopping E

% Initial E/V/I calculation
[V,I,E,T_v,dT]=RECAL_i(VI,C_diag,s);
% Initial Probability setting
[Probability,total_rate]=Update(E,T,OV_list,T_hk);

%% MMC Experiments

criteria=1;
counter=0;
time=0;
mini_batch=20;
total_time=0;
converge_count=0;
[E_value,E_position]=max(E(:));
Check_Volum=[];
Entropy_list=[];
Result_list=[0,0,I,E_value,E_position];
I_list=[I];
NC_list=cell(1,4);

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

            temp_value1=total_rate(:,pick_y+2)-seed(i,2);
            temp_value1(temp_value1<0)=inf;
            [~,pick_x]=min(temp_value1);
            process_Probability=Probability{pick_x,pick_y};
            if ~isempty(process_Probability)
                fix=0;
            else
                seed(i,:)=rand(1,3);
            end
        end
        temp_value2=process_Probability(:,1)-seed(i,3);
        temp_value2(temp_value2<0)=inf;            
        [~,ind]=min(temp_value2);
        
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
        beta=1/(log(sum(total_rate(:))));
        delt_time(i)=(log(1/(seed(i,3)^3)))*(beta)*(n1*n2*aphla_C^2/(n_w*n_l*n_f))^(1/log(sum(total_rate(:))/(n1*n2*(n3-2))))*(VI^(v_factor_1)*exp(v_factor_2)/10);
        
        batch_time=batch_time+delt_time(i);
        
    end
    %% Decide Entropy for diverse events
    
    [V,temp_I,temp_E,T_v,T,OV_list,temp_Probability,temp_total_rate,criteria_index,fb_list,check_volum_f,total_entropy,prob]=MMC_decide(defect_configure,mini_batch,VI,V,T_v,T,OV_list,T_hk);
    Entropy_list=[Entropy_list;total_entropy,prob,sum(dT(:))-sum(T(:))];
    
    if ~isempty(temp_I)
        I=temp_I;
        E=temp_E;
        Probability=temp_Probability;
        total_rate=temp_total_rate;
        total_time=total_time+batch_time;
        if criteria_index==1
            criteria=1;
            Check_Volum=check_volum_f;
        end
%         converge_count=0; % Set 0 for not converge
        counter=counter+1;
        [E_value,E_position]=max(E(:));
        Result_list=[Result_list;total_time,sum(M(:)),I,E_value,E_position];
        I_list=[I_list,I];
        converge_count=0;
    else
        [Probability,total_rate]=Update(E,T,OV_list,T_hk);
        converge_count=converge_count+1;
        if converge_count>=10 % continous 5 minibatch as converge criteria
            criteria=0;
            fprintf('\n Convergence Approaching')
        end
        
    end
    
    if counter == 20
        NC_list{1} = M;
    elseif counter==40
        NC_list{2}= M;
    elseif counter==60
        NC_list{3}=M;
    elseif counter==80
        NC_list{4}=M;
    end
    
    fprintf('\n Counter: %d; Total Time: %e; I: %e;',counter,total_time,I)
%     if I>Result_list(1,3)*scale*3
    if counter==0
        continue
    else
        if (I/I_list(end-1)>2)&&(length(I_list)>5)
            criteria=0;
            fprintf('\n Breakdown')
        end
    end
end

% Result_list(:,3)=Result_list(:,3)/Result_list(end-1,3)*I_factor;
% Result_list(:,1)=Result_list(:,1)/Result_list(end,1)*t_factor;
% Result_list(:,2)=Result_list(:,2)/Result_list(end,2)*t_factor;
% 
% k=D{chose_x,chose_y}(:,2);
% k(end)=Result_list(end,3);
% figure
% hold on
% plot(D{chose_x,chose_y}(:,1),k)
% plot(Result_list(:,1),Result_list(:,3))
% hold off


% Result_list(:,1)=Result_list(:,1)/Result_list(end,1)*t_factor;
% Result_list(:,2)=Result_list(:,2)/Result_list(end,2)*t_factor;



I_s=Result_list(:,3)*(n_w*n_l*n_f)/(n1*n2*aphla_C^2);
Result_list(:,3)=I_s;


s=Result_list(:,3)/Result_list(3,3)*Initial_Ipre*Result_list(2,1)^(1.5*slope);


figure
hold on
% plot(D{chose_x,chose_y}(:,1),D{chose_x,chose_y}(:,2))
plot(Result_list(:,1),s)
hold off

