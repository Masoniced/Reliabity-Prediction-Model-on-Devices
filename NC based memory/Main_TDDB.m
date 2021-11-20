%% Main_TDDB

load('C:\D\Program Files\Project\Matlab Project\Degradation\Overall\Initializing Para\Initial_125_70_0.4.mat')
% Data=load('C:\Users\Meisen\Documents\MATLAB\projeects\Clustering data processing\Fourier\current.mat');
% D=Data.res;
% chose_x=2;
% chose_y=5;
% scale=D{chose_x,chose_y}(end-1,2)/D{chose_x,chose_y}(1,2);
% I_factor=D{chose_x,chose_y}(end-2,2);
% t_factor=D{chose_x,chose_y}(end,1);

VI=1;

n_w=7e-8;
n_l=1e-7;
n_f=22;
t_o=3.8e-9;

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


number=1;
final_list={};

while number<=200

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
                criteria=0;
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

    %     fprintf('\n Counter: %d; Total Time: %e; I: %e',counter,total_time,I)
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
    
    temp_current=Result_list(:,3)/Result_list(3,3)*Initial_Ipre*Result_list(3,1)^(0.8*slope);
    Result_list(:,3)=temp_current;
    final_list{number,1}=[Result_list(end,1),Result_list(end,2),Result_list(end,3),Result_list(end,4)];
    final_list{number,2}=Result_list;

    fprintf('\n Number: %d; Total Time: %e; I: %e',number,final_list{number,1}(1),final_list{number,1}(3))
    load('C:\D\Program Files\MATLAB\projeects\Degradation\Overall\Initializing Para\Initial_125_70_0.4.mat')
    
    VI=1;
    Initial_Ipre=normrnd(Initial_Ipre_mean,viaration_Ipre);

    Ea_s=0.7*e; % J
    Er_h=1.0*e; % HK Recombination E
    Er_s=1.95*e; % IL Recombination E
    Eh_h=1.1*e; % HK Hopping E
    Eh_s=1.9*e; % IL Hopping E
    number=number+1;
end


