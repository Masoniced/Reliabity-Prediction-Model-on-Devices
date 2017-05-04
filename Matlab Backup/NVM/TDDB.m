%% MCMC percolation path sequence
syms samp
global M Lable_mat label Label_record Rate n1 n2 R
global Ks e Ea_s Ea_h p_s p_h kh kh_gb ks_v ks_p K S_v M_s
%% Basic parameter setting
n1=100;n2=8; % assume each cell as 0.8A, thus HK 4.8nm, IL 1.6nm; overall diameter 80*80*6.4
M=zeros(n1,2,n2); %occupation matrix
Rate=zeros(n1,2,n2); %rate matrix
%% Cluster checking parameter
Label_record={};
Lable_mat=zeros(size(M));
label=0;
%% Specific parameter setting
% Model constant
M_s=8e-10; % m
% Typical vibration frequency
Ks=1e13; % Hz
% Electron charge
e=1.6e-19; % C
% Free activation of IL(SiOx)
Ea_s=2.43*e; % J
% Free activation of HK(HfO2)
Ea_h=4.4*e; % J
% DM of IL
p_s=4.875*e*1e-10; % e*`A
% DM of HK
p_h=10.2*e*1e-10; % e*`A
% Permitt of HK
kh=25;
% Permitt of HK GB
kh_gb=26.5;
% Permitt of virgin IL
ks_v=6;
% Permitt of Post BD IL
ks_p=8.5;
% Boltzmann constant
K=1.38e-23;
% Stressing voltage
S_v=2; % V
%% KMC
criteria=0;
n_gb=2; % number of random GB path
PIT_p=0.001; % percentage of process introduced failure
[R,compare_volum]=RPG_r(criteria,PIT_p);
U=load('C:\D\Program Files\MATLAB\projeects\TDDB\Data\matlab.mat');
R=U.R;
Rate_vector_HK=zeros(n1*2*(n2-2),2);
Rate_vector_IL=zeros(n1*2*2,2);
total_rate=0;
es_kh=kh;
es_ks=ks_v;
E_field_hk=S_v*(1/((es_kh/6)*(1/es_ks)+1))/(M_s*6);
E_field_hk_ind=S_v*(1/((es_kh/6)*(1/ks_p)+1))/(M_s*6);
E_field_il=S_v*(1/((6/es_kh)*(es_ks/1)+1))/(M_s*2);
sample_rate_hk1=Ks*exp(-(Ea_h-p_h*(kh+2)/3*E_field_hk_ind)/(K*300));
sample_rate_hk2=Ks*exp(-(Ea_h-p_h*(kh+2)/3*E_field_hk_ind)/(K*300))*1.5;
sample_rate_hk3=Ks*exp(-(Ea_h-p_h*(kh_gb+2)/3*E_field_hk_ind)/(K*300));
sample_rate_hk4=Ks*exp(-(Ea_h-p_h*(kh_gb+2)/3*E_field_hk_ind)/(K*300))*1.5;
sample_rate_il1=Ks*exp(-(Ea_s-p_s*(ks_v+2)/3*E_field_il)/(K*300));
sample_rate_il2=Ks*exp(-(Ea_s-p_s*(ks_p+2)/3*E_field_il)/(K*300));
Rate_vector_IL_pre=0;
Rate_vector_HK_pre=0;
for i=1:1:n1*2*n2
    [a,b,c]=ind2sub([n1,2,n2],i);
    if c==n2||c==n2-1
        if R(a,b,c)==1
            Rate(a,b,c)=sample_rate_il2;
        else
            Rate(a,b,c)=sample_rate_il1;
        end
        p_p=sub2ind([n1,2,2],a,b,c-n2+2);
        Rate_vector_IL(p_p,1)=Rate(a,b,c)+Rate_vector_IL_pre;
        Rate_vector_IL(p_p,2)=p_p;
        Rate_vector_IL_pre=Rate_vector_IL(p_p,1);
    else
        if R(a,b,c)==1
            Rate(a,b,c)=sample_rate_hk3;
        else
            Rate(a,b,c)=sample_rate_hk1;
        end
        p_p=sub2ind([n1,2,n2-2],a,b,c);
        Rate_vector_HK(p_p,1)=Rate(a,b,c)+Rate_vector_HK_pre;
        Rate_vector_HK(p_p,2)=p_p;
        Rate_vector_HK_pre=Rate_vector_HK(p_p,1);
    end
end
Rate_vector_IL(:,1)=Rate_vector_IL(:,1)/Rate_vector_IL_pre;
Rate_vector_HK(:,1)=Rate_vector_HK(:,1)/Rate_vector_HK_pre;
%% Experiment
counter=0;
time=0;
total_hk_rate=Rate_vector_HK_pre;
total_il_rate=Rate_vector_IL_pre;
change_position=zeros(1,3);
so=0;
check_volum=[];
check_value=Rate_vector_HK(1,1);
il_count=0;
il_time=0;
S=zeros(n1,2,n2);
fb=0;
criteria_1=0;
check1=0;
check2=0;
while criteria==0
    rng('shuffle')
    f_hk=rand(1,1);
    f_il=rand(1,1);
    counter=counter+1;
    if counter==n1*n1*2
        break
    end
    % find the position
    [~,hk_p]=min(abs(Rate_vector_HK(:,1)-f_hk));
    [~,il_p]=min(abs(Rate_vector_IL(:,1)-f_il));
    [p_hk(1),p_hk(2),p_hk(3)]=ind2sub([n1,2,n2-2],Rate_vector_HK(hk_p,2));
    [p_il(1),p_il(2),p_il(3)]=ind2sub([n1,2,2],Rate_vector_IL(il_p,2));
    p_il(3)=p_il(3)+n2-2;
    % p_hk=[Rate_vector_HK(hk_p,2),Rate_vector_HK(hk_p,3),Rate_vector_HK(hk_p,4)];
    % p_il=[Rate_vector_IL(il_p,2),Rate_vector_IL(il_p,3),Rate_vector_IL(il_p,4)];
    % time increase
    time=time+log(1/f_hk)*(1/total_hk_rate)*(n1*2*(n2-2)/(n1*2*(n2-2)-counter))+log(1/f_il)*(1/total_il_rate)*(n1*2*(2)/(n1*2*(2)-counter));
    % update occupy matrix
    M(p_hk(1),p_hk(2),p_hk(3))=1;
    M(p_il(1),p_il(2),p_il(3))=1;
    % Checking cluster and criteria
    [check1,~,decide_label1]=Cluster(p_hk);
    if check1==0
        [check2,fb,decide_label2]=Cluster(p_il);
    end
    if (check1==1)&&(check2~=1)
        criteria_1=1;
        decide_label=decide_label1;
    elseif (check2==1)&&(check1~=1)
        criteria_1=1;
        decide_label=decide_label2;
    elseif (check1==1)&&(check2==1)
        criteria_1=1;
        decide_label=decide_label2;
    else
        criteria_1=0;
    end
    % update of dynamic values (total_hk_rate,total_il_rate,Rate)
    if fb==1
        il_count=il_count+1;
        if p_il(3)==n2-1
            change_position(1)=p_il(1);
            change_position(2)=p_il(2);
            change_position(3)=p_il(3)-1;
        else
            change_position(1)=p_il(1);
            change_position(2)=p_il(2);
            change_position(3)=p_il(3)-2;
        end
        if il_count==1
            S=M;
            il_time=time;
        end
    else
        change_position(1)=0;
        change_position(2)=0;
        change_position(3)=0;
    end
    temp_hk_rate=total_hk_rate;
    temp_il_rate=total_il_rate;
    if (fb==1)&&(M(change_position(1),change_position(2),change_position(3))==0)&&(p_hk(1)~=change_position(1)||p_hk(2)~=change_position(2)||p_hk(3)~=change_position(3))
        if R(change_position(1),change_position(2),change_position(3))==1
            total_hk_rate=total_hk_rate+sample_rate_hk4-sample_rate_hk3-Rate(p_hk(1),p_hk(2),p_hk(3));
            Rate(change_position(1),change_position(2),change_position(3))=sample_rate_hk4;
        else
            total_hk_rate=total_hk_rate+sample_rate_hk2-sample_rate_hk1-Rate(p_hk(1),p_hk(2),p_hk(3));
            Rate(change_position(1),change_position(2),change_position(3))=sample_rate_hk2;
        end
        total_il_rate=total_il_rate-Rate(p_il(1),p_il(2),p_il(3));
        so=1;
    else
        total_hk_rate=total_hk_rate-Rate(p_hk(1),p_hk(2),p_hk(3));
        total_il_rate=total_il_rate-Rate(p_il(1),p_il(2),p_il(3));
        so=0;
    end
    % update rate vectors
    Rate_vector_HK(hk_p,:)=[];
    if M(p_hk(1),p_hk(2),p_hk(3))==0
        fprintf('update rate vectors error')
        break
    end
    Rate_vector_IL(il_p,:)=[];
    if M(p_il(1),p_il(2),p_il(3))==0
        fprintf('update rate vectors error')
        break
    end    
    % For HK
    switch so
        case 1
            pre_change_ind=sub2ind([n1,n1,n2-2],change_position(1),change_position(2),change_position(3));
            change_ind=find(Rate_vector_HK(:,2)==pre_change_ind);
            if length(change_ind)~=1
                break
            end
            if hk_p==(n1*n1*(n2-2)-counter+1)
                if change_ind==1
                    Rate_vector_HK(1:1:end,1)=Rate_vector_HK(1:1:end,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                elseif change_ind==(n1*n1*(n2-2)-counter)
                    Rate_vector_HK(1:1:(change_ind-1),1)=Rate_vector_HK(1:1:(change_ind-1),1)*(temp_hk_rate)/(total_hk_rate);
                    Rate_vector_HK(change_ind,1)=Rate_vector_HK(change_ind,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                else
                    Rate_vector_HK(1:1:(change_ind-1),1)=Rate_vector_HK(1:1:(change_ind-1),1)*(temp_hk_rate)/(total_hk_rate);
                    Rate_vector_HK(change_ind:1:end,1)=Rate_vector_HK(change_ind:1:end,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                end
%                if (total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))<0
%                    fprintf('Update vector error 1')
%                    return
%                end
            elseif hk_p==1
                Rate_vector_HK(1:1:end,1)=Rate_vector_HK(1:1:end,1)-Rate(p_hk(1),p_hk(2),p_hk(3))/(temp_hk_rate);
                if change_ind==1
                    Rate_vector_HK(1:1:end,1)=Rate_vector_HK(1:1:end,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                elseif change_ind==(n1*n1*(n2-2)-counter)
                    Rate_vector_HK(1:1:(change_ind-1),1)=Rate_vector_HK(1:1:(change_ind-1),1)*(temp_hk_rate)/(total_hk_rate);
                    Rate_vector_HK(change_ind,1)=Rate_vector_HK(change_ind,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                else
                    Rate_vector_HK(1:1:(change_ind-1),1)=Rate_vector_HK(1:1:(change_ind-1),1)*(temp_hk_rate)/(total_hk_rate);
                    Rate_vector_HK(change_ind:1:end,1)=Rate_vector_HK(change_ind:1:end,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                end
%                if (total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))<0
%                    fprintf('Update vector error 2')
%                    return
%                end
            else
                Rate_vector_HK(hk_p:1:end,1)=Rate_vector_HK(hk_p:1:end,1)-Rate(p_hk(1),p_hk(2),p_hk(3))/(temp_hk_rate);
                if change_ind==1
                    Rate_vector_HK(1:1:end,1)=Rate_vector_HK(1:1:end,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                elseif change_ind==(n1*n1*(n2-2)-counter)
                    Rate_vector_HK(1:1:(change_ind-1),1)=Rate_vector_HK(1:1:(change_ind-1),1)*(temp_hk_rate)/(total_hk_rate);
                    Rate_vector_HK(change_ind,1)=Rate_vector_HK(change_ind,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                else
                    Rate_vector_HK(1:1:(change_ind-1),1)=Rate_vector_HK(1:1:(change_ind-1),1)*(temp_hk_rate)/(total_hk_rate);
                    Rate_vector_HK(change_ind:1:end,1)=Rate_vector_HK(change_ind:1:end,1)*(temp_hk_rate)/(total_hk_rate)+(total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))/(total_hk_rate);
                end
%                if (total_hk_rate-temp_hk_rate+Rate(p_hk(1),p_hk(2),p_hk(3)))<0
%                    fprintf('Update vector error 3')
%                    return
%                end                
            end
        case 0
            if hk_p==(n1*n1*(n2-2)-counter+1)
                Rate_vector_HK(1:1:end,1)=Rate_vector_HK(1:1:end,1)*(temp_hk_rate)/(total_hk_rate);
            elseif hk_p==1
                Rate_vector_HK(1:1:end,1)=Rate_vector_HK(1:1:end,1)*(temp_hk_rate)/(total_hk_rate)-Rate(p_hk(1),p_hk(2),p_hk(3))/(total_hk_rate);
            else
                Rate_vector_HK(1:1:hk_p-1,1)=Rate_vector_HK(1:1:hk_p-1,1)*(temp_hk_rate)/(total_hk_rate);
                Rate_vector_HK(hk_p:1:end,1)=Rate_vector_HK(hk_p:1:end,1)*(temp_hk_rate)/(total_hk_rate)-Rate(p_hk(1),p_hk(2),p_hk(3))/(total_hk_rate);
            end
        otherwise
            fprintf('case error')
            fprintf('%d\n',s,p_hk,p_il)
    end
    % For IL
    if il_p==n1*n1*2-counter+1
        Rate_vector_IL(1:1:il_p-1,1)=Rate_vector_IL(1:1:il_p-1,1)*(temp_il_rate)/(total_il_rate);
    elseif il_p==1
        Rate_vector_IL(il_p:1:end,1)=Rate_vector_IL(il_p:1:end,1)*(temp_il_rate)/(total_il_rate)-Rate(p_il(1),p_il(2),p_il(3))/(total_il_rate);
    else
        Rate_vector_IL(il_p:1:end,1)=Rate_vector_IL(il_p:1:end,1)*(temp_il_rate)/(total_il_rate)-Rate(p_il(1),p_il(2),p_il(3))/(total_il_rate);
        Rate_vector_IL(1:1:il_p-1,1)=Rate_vector_IL(1:1:il_p-1,1)*(temp_il_rate)/(total_il_rate);
    end
    if Rate_vector_HK(1,1)/check_value>1e2
        fprintf('IL case error')
        break
    else
        check_value=Rate_vector_HK(1,1);
    end
    % Output for real-time checking
    if mod(counter,100)==0
        fprintf('%d\n',counter,time)
    end
    if criteria_1==1
        vector=[Label_record{:,2}];
        check_vector=find(vector==decide_label);
        w=length(check_vector);
        for i=1:w
            check_volum=[check_volum;Label_record{check_vector(i),1}];
        end
        uu=length(check_volum(:,1));
        point=[sum(check_volum(:,1)/uu),sum(check_volum(:,2)/uu),0];
        point(1)=(point(1)-0.5-n1/2)*0.8;
        point(2)=(point(2)-0.5-n1/2)*0.8;
    end
    if criteria_1==1
        criteria=1;
    end
end
weibull=2*counter; 
num1=length(check_volum);
s_L=zeros(num1,3);
u=1./check_volum;
check_number=length(find(compare_volum*u'==3));
    
    