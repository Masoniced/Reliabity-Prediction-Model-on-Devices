function [V,I,E,T_v,T,OV_list,Probability,total_rate,criteria_index,fb_list,check_volum_f,total_entropy,prob]=MMC_decide(defect_configure,num_conf,VI,V,T_v,T,OV_list,T_hk,scale)
% structure of defect configuration : cell(num_cof*3) 
% 1:event type;
% 2:position information;
% structure need to be modified:tunneling_lable,V,T_v,T
%
global location_conduct tunneling_lable Thermal_contribution Location_lable ovcp_IL ovcp_HK Entro_HK Entro_IL n1 n2 n3 M K

O_concentration=1e-6;
HK_Add_volum=[];
IL_Add_volum=[];
HK_Remove_volum=[];
IL_Remove_volum=[];
C_add_volum=[];
C_remove_volum=[];
Thermal_add_volum=[];
Thermal_remove_volum=[];

for i=1:num_conf
    s=defect_configure{i,1};
    switch s
        case 1
            if defect_configure{i,3}=='H'
                HK_Add_volum=[HK_Add_volum,defect_configure{i,2}];
            else
                IL_Add_volum=[IL_Add_volum,defect_configure{i,2}];
            end
            C_add_volum=[C_add_volum,location_conduct{defect_configure{i,2}}];
            Thermal_add_volum=[Thermal_add_volum,Location_lable(defect_configure{i,2},:)];
        case 2
            if defect_configure{i,3}=='H'
                HK_Remove_volum=[HK_Remove_volum,defect_configure{i,2}];
            else
                IL_Remove_volum=[IL_Remove_volum,defect_configure{i,2}];
            end
            C_remove_volum=[C_remove_volum,location_conduct{defect_configure{i,2}}];
            Thermal_remove_volum=[Thermal_remove_volum,Location_lable(defect_configure{i,2},:)];
        case 3
            if defect_configure{i,3}(2)=='H'
                HK_Add_volum=[HK_Add_volum,defect_configure{i,2}(2)];
            else
                IL_Add_volum=[IL_Add_volum,defect_configure{i,2}(2)];
            end
            if defect_configure{i,3}(1)=='H'
                HK_Remove_volum=[HK_Remove_volum,defect_configure{i,2}(1)];
            else
                IL_Remove_volum=[IL_Remove_volum,defect_configure{i,2}(1)];
            end            
            C_add_volum=[C_add_volum,location_conduct{defect_configure{i,2}(2)}];
            C_remove_volum=[C_remove_volum,location_conduct{defect_configure{i,2}(1)}];
            
            Thermal_add_volum=[Thermal_add_volum,Location_lable(defect_configure{i,2}(2),:)];
            Thermal_remove_volum=[Thermal_remove_volum,Location_lable(defect_configure{i,2}(1),:)];
    end
end
%% Update of candidates
L_HK_A=length(HK_Add_volum);
L_HK_R=length(HK_Remove_volum);
L_IL_A=length(IL_Add_volum);
L_IL_R=length(IL_Remove_volum);

if ~isempty(Thermal_add_volum)
    Thermal_add_volum=sort(Thermal_add_volum);
    Thermal_add_volum=Thermal_add_volum([true,diff(Thermal_add_volum)>0]);
end

if ~isempty(Thermal_remove_volum)
Thermal_remove_volum=sort(Thermal_remove_volum);
Thermal_remove_volum=Thermal_remove_volum([true,diff(Thermal_remove_volum)>0]);
end

previous_defect_number=length(OV_list{1})+length(OV_list{2})+1;

Thermal_contribution(Thermal_contribution>1)=previous_defect_number+L_HK_A+L_IL_A;
Thermal_contribution(Thermal_add_volum)=previous_defect_number+L_HK_A+L_IL_A;
Thermal_contribution(Thermal_remove_volum)=1;


C_add_volum=sort(C_add_volum);
C_remove_volum=sort((C_remove_volum));
double_add=C_add_volum(~[true,diff(C_add_volum)>0]);
double_remove=C_remove_volum(~[true,diff(C_remove_volum)>0]);

tunneling_lable(C_add_volum)=tunneling_lable(C_add_volum)+1;
tunneling_lable(double_add)=tunneling_lable(double_add)+1;
tunneling_lable(C_remove_volum)=tunneling_lable(C_remove_volum)-1;
tunneling_lable(double_remove)=tunneling_lable(double_remove)-1;


if nargin==9
    [C_diag,s]=update_C(V,T_v,VI,scale);
else
    [C_diag,s]=update_C(V,T_v,VI);
end

[temp_V,temp_I,temp_E,temp_T_v,dT]=RECAL_i(VI,C_diag,s);
% if temp_I>=1e-6
%     fprintf('over current')
%     return
% elseif temp_I<0
%     fprintf('negative current')
%     return
% end

sum_HKdT=sum(dT([HK_Add_volum,HK_Remove_volum]));
sum_HKT=sum(T([HK_Add_volum,HK_Remove_volum]));
sum_ILdT=sum(dT([IL_Add_volum,IL_Remove_volum]));
sum_ILT=sum(T([IL_Add_volum,IL_Remove_volum]));
criteria_index=0;
check_volum_f=[];
fb_list=zeros(L_HK_A+L_IL_A,1);
%% Calculate maximum entropy
% total_entropy=ovcp_HK*O_concentration*(L_HK_A-L_HK_R)....
%     +ovcp_IL*O_concentration*(L_IL_A-L_IL_R)....
%     -Entro_HK*(sum_HKdT-sum_HKT)....
%     -Entro_IL*(sum_ILdT-sum_ILT);
total_entropy=ovcp_HK*O_concentration*(L_HK_A-L_HK_R)....
    +ovcp_IL*O_concentration*(L_IL_A-L_IL_R)....
    -Entro_HK*(sum(dT(1:n1*n2*T_hk))-sum(T(1:n1*n2*T_hk)))....
    -Entro_IL*(sum(dT(n1*n2*T_hk+1:end))-sum(T(n1*n2*T_hk+1:end)));
if total_entropy<=0
    % Basic Update
    V=temp_V;
    I=temp_I;
    E=temp_E;
    T_v=temp_T_v;
    T=dT;   
    % HK OV_list Update
    ind_HK=ismembc([OV_list{1}],sort(HK_Remove_volum));
    OV_list{1}=OV_list{1}(~ind_HK);
    OV_list{1}=[OV_list{1},HK_Add_volum];
    % IL OV_list Update
    IL_Remove_volum=IL_Remove_volum-n1*n2*T_hk;
    IL_Add_volum=IL_Add_volum-n1*n2*T_hk;
    ind_IL=ismembc([OV_list{2}],sort(IL_Remove_volum));
    OV_list{2}=OV_list{2}(~ind_IL);
    OV_list{2}=[OV_list{2},IL_Add_volum];
    % Probability Vector Update
    [Probability,total_rate]=Update(E,T,OV_list,T_hk);
    % Clustering Update
    for i=1:L_HK_R
        f=HK_Remove_volum(i);
        p=zeros(1,3);
        p(1)=rem(f,n1);
        if p(1)==0
            p(1)=n1;
        end
        p(2)=rem((f-p(1))/n1,n2)+1;
        p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
        M(f)=0;
        [~,~,~,~]=Cluster_AR(p,1,n1,n2,n3);
    end
    for i=1:L_IL_R
        f=IL_Remove_volum(i)+n1*n2*T_hk;
        p=zeros(1,3);
        p(1)=rem(f,n1);
        if p(1)==0
            p(1)=n1;
        end        
        p(2)=rem((f-p(1))/n1,n2)+1;
        p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
        M(f)=0;
        [~,~,~,~]=Cluster_AR(p,1,n1,n2,n3);
    end       
    for i=1:L_HK_A
        f=HK_Add_volum(i);
        p=zeros(1,3);
        p(1)=rem(f,n1);
        if p(1)==0
            p(1)=n1;
        end        
        p(2)=rem((f-p(1))/n1,n2)+1;
        p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
        M(f)=1;
        [criteria,fb,~,check_volum_find]=Cluster_AR(p,0,n1,n2,n3);
        if criteria==1
            criteria_index=1;
            check_volum_f=check_volum_find;
        end
        fb_list(i)=fb;
    end             
    for i=1:L_IL_A
        f=IL_Add_volum(i)+n1*n2*T_hk;
        p=zeros(1,3);
        p(1)=rem(f,n1);
        if p(1)==0
            p(1)=n1;
        end        
        p(2)=rem((f-p(1))/n1,n2)+1;
        p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
        M(f)=1;
        [criteria,fb,~,check_volum_find]=Cluster_AR(p,0,n1,n2,n3);
        if criteria==1
            criteria_index=1;
            check_volum_f=check_volum_find;
        end
        fb_list(i+L_HK_A)=fb;    
    end
    prob=0;
elseif total_entropy>0
    seed=prod(rand(1,10));
    prob=exp(-total_entropy*(L_HK_A+L_HK_R+L_IL_A+L_IL_R)/(K*(sum_HKdT+sum_ILdT)));
    if seed<1-prob
        % Basic Update
        V=temp_V;
        I=temp_I;
        E=temp_E;
        T_v=temp_T_v;
        T=dT;        
        % HK OV_list Update
        ind_HK=ismembc([OV_list{1}],sort(HK_Remove_volum));
        OV_list{1}=OV_list{1}(~ind_HK);
        OV_list{1}=[OV_list{1},HK_Add_volum];
        % IL OV_list Update
        IL_Remove_volum=IL_Remove_volum-n1*n2*T_hk;
        IL_Add_volum=IL_Add_volum-n1*n2*T_hk;
        ind_IL=ismembc([OV_list{2}],sort(IL_Remove_volum));
        OV_list{2}=OV_list{2}(~ind_IL);
        OV_list{2}=[OV_list{2},IL_Add_volum];
        % Probability Vector Update
        [Probability,total_rate]=Update(E,T,OV_list,T_hk);
        % Clustering Update
        for i=1:L_HK_R
            f=HK_Remove_volum(i);
            p=zeros(1,3);
            p(1)=rem(f,n1);
            if p(1)==0
                p(1)=n1;
            end
            p(2)=rem((f-p(1))/n1,n2)+1;
            p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
            M(f)=0;
            [~,~,~,~]=Cluster_AR(p,1,n1,n2,n3);
        end
        for i=1:L_IL_R
            f=IL_Remove_volum(i)+n1*n2*T_hk;
            p=zeros(1,3);
            p(1)=rem(f,n1);
            if p(1)==0
                p(1)=n1;
            end            
            p(2)=rem((f-p(1))/n1,n2)+1;
            p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
            M(f)=0;
            [~,~,~,~]=Cluster_AR(p,1,n1,n2,n3);
        end       
        for i=1:L_HK_A
            f=HK_Add_volum(i);
            p=zeros(1,3);
            p(1)=rem(f,n1);
            if p(1)==0
                p(1)=n1;
            end            
            p(2)=rem((f-p(1))/n1,n2)+1;
            p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
            M(f)=1;
            [criteria,fb,~,check_volum_find]=Cluster_AR(p,0,n1,n2,n3);
            if criteria==1
                criteria_index=1;
                check_volum_f=check_volum_find;
            end
            fb_list(i)=fb;            
        end             
        for i=1:L_IL_A
            f=IL_Add_volum(i)+n1*n2*T_hk;
            p=zeros(1,3);
            p(1)=rem(f,n1);
            if p(1)==0
                p(1)=n1;
            end            
            p(2)=rem((f-p(1))/n1,n2)+1;
            p(3)=(f-p(1)-n1*(p(2)-1))/(n1*n2)+1;
            M(f)=1;
            [criteria,fb,~,check_volum_find]=Cluster_AR(p,0,n1,n2,n3);
            if criteria==1
                criteria_index=1;
                check_volum_f=check_volum_find;
            end
            fb_list(i+L_HK_A)=fb;              
        end  
    else
        tunneling_lable(C_add_volum)=tunneling_lable(C_add_volum)-1;
        tunneling_lable(double_add)=tunneling_lable(double_add)-1;
        tunneling_lable(C_remove_volum)=tunneling_lable(C_remove_volum)+1;
        tunneling_lable(double_remove)=tunneling_lable(double_remove)+1;
        
        Thermal_contribution(Thermal_add_volum)=1;
        Thermal_contribution(Thermal_contribution>1)=previous_defect_number;
        Thermal_contribution(Thermal_remove_volum)=previous_defect_number;
        
        I=[];
        E=[];
        Probability=[];
        total_rate=[];
        criteria_list=[];
        fb_list=[];
    end
end
% if length([Label_record{:,2}])~=sum(M(:))
%     fprintf('\n Counting Error')
%     return
% end
end