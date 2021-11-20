
%% Probability Vector Update
function [Probability,total_rate]=Update(E,T,OV_list,T_hk)
global n1 n2 n3 K Ks D_K % Basic 
global Ea_h Ea_s p_h p_s % Generation
global Er_h Er_s % Recombination
global Eh_h Eh_s % Hopping
global Initial_T
%tic
%% Output Demonstration
%     Probability: 3*2                               
%  HK_G_vector  IL_G_vector       
%  HK_R_vector  IL_R_vector
%  HK_H_vector  IL_H_vector

%                            total_rate: 3*5
%  HK_G_total   HK_G_total   HK_G_proportion   IL_G_proportion   Total_HK
%  HK_R_total   HK_R_total   HK_R_proportion   IL_R_proportion   Total_HK
%  HK_H_total   HK_H_total   HK_H_proportion   IL_H_proportion   Total_HK
%% Basic Setting

s=2;
s1=2;
pf=2.5;

hkl=length(OV_list{1});
ill=length(OV_list{2});
E_HK=E(:,:,1:T_hk);
E_HK=E_HK/2;
E_IL=E(:,:,T_hk+1:end);
T_HK=T(:,:,1:T_hk);
T_IL=T(:,:,T_hk+1:end);
KH=D_K(:,:,1:T_hk);
KS=D_K(:,:,T_hk+1:end);

% C_OX_HK=C_OX(:,:,1:T_hk);
% C_OX_IL=C_OX(:,:,T_hk+1:end);


factor_r_HK=Ks*hkl^8/(n1*n2*T_hk);
factor_r_IL=Ks*(ill^3/(n1*n2*(n3-2-T_hk)));
factor_h_HK=hkl^8/log(n1*n2*T_hk);
factor_h_IL=ill^3/log(n1*n2*(n3-2-T_hk));

total_rate=zeros(3,5);
c1=0;
c2=0;

%% Generation Update

HK_G_vector=zeros(n1*n2*T_hk,3);
HK_G_vector(:,2)=1:n1*n2*T_hk;
IL_G_vector=zeros(n1*n2*(n3-2-T_hk),3);
IL_G_vector(:,2)=1:n1*n2*(n3-2-T_hk);


HK_G_probability=Ks*exp(-Ea_h./(K*T_HK)).*exp(p_h*(KH+2)/3.*E_HK./(K*T_HK));
IL_G_probability=Ks*exp(-Ea_s./(K*T_IL)).*exp(p_s*(KS+2)/3.*E_IL./(K*T_IL));

HK_G_probability=reshape(HK_G_probability,[n1*n2*T_hk,1]);
IL_G_probability=reshape(IL_G_probability,[n1*n2*(n3-2-T_hk),1]);

HK_G_vector(:,3)=HK_G_probability;
IL_G_vector(:,3)=IL_G_probability;

if ~isempty(OV_list{1})
    HK_G_probability(OV_list{1})=[];
    HK_G_vector(OV_list{1},:)=[];
end
if ~isempty(OV_list{2})
    IL_G_probability(OV_list{2})=[];
    IL_G_vector(OV_list{2},:)=[];
end
if isempty(OV_list{1})
    c1=1;
end
if isempty(OV_list{2})
    c2=1;
end

total_rate(1,1)=sum(HK_G_probability);
total_rate(1,2)=sum(IL_G_probability);

HK_G_probability=HK_G_probability/(min(HK_G_probability))+s;
IL_G_probability=IL_G_probability/(min(IL_G_probability))+s;

HK_G_probability=cumsum(log(HK_G_probability).^pf);
IL_G_probability=cumsum(log(IL_G_probability).^pf);

HK_G_vector(:,1)=HK_G_probability/HK_G_probability(end);
IL_G_vector(:,1)=IL_G_probability/IL_G_probability(end);

%% Recombination Update && Hopping Update
if c1==1&&c2==1
    
    HK_H_vector=[];
    IL_H_vector=[];  
    total_rate(2,1)=0;
    total_rate(2,2)=0;
    total_rate(3,1)=0;
    total_rate(3,2)=0;
    HK_R_vector=[];
    IL_R_vector=[];
    
elseif c1==1&&c2==0
    
    IL_R_vector=zeros(ill,3);
    IL_R_vector(:,2)=OV_list{2};
    
    IL_R_probability=factor_r_IL*exp((T_IL(OV_list{2})-Initial_T)/30+3).*exp(-Er_s./(K*T_IL(OV_list{2})));
    IL_R_probability=IL_R_probability';
    IL_R_vector(:,3)=IL_R_probability;
    
    total_rate(2,1)=0;
    total_rate(2,2)=sum(IL_R_probability);
    HK_R_vector=[];
    
    IL_R_probability=IL_R_probability/(min(IL_R_probability))+s;
    
    IL_R_probability=cumsum(log(IL_R_probability).^pf);
    
    IL_R_vector(:,1)=IL_R_probability/IL_R_probability(end);

    % Recombination
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hopping
    
    IL_H_vector=zeros(ill,3);
    IL_H_vector(:,2)=OV_list{2}; 
    
    IL_H_probability=Ks*factor_h_IL*exp(-Eh_s./(K*T_IL(OV_list{2}))).*exp(p_s*(KS(OV_list{2})+2)/3.*E_IL(OV_list{2})./(K*T_IL(OV_list{2})));
    IL_H_probability=IL_H_probability';
    IL_H_vector(:,3)=IL_H_probability;
    

    total_rate(3,1)=0;
    total_rate(3,2)=sum(IL_H_probability);   
    HK_H_vector=[];
    
    IL_H_probability=IL_H_probability/(min(IL_H_probability))+s;
    
    IL_H_probability=cumsum(log(IL_H_probability).^pf);    
    
    IL_H_vector(:,1)=IL_H_probability/IL_H_probability(end);  
    
elseif c2==1&&c1==0
    
    HK_R_vector=zeros(hkl,3);
    HK_R_vector(:,2)=OV_list{1};
    
    HK_R_probability=factor_r_HK*exp((T_HK(OV_list{1})-Initial_T)/10+3).*exp(-Er_h./(K*T_HK(OV_list{1})));
    HK_R_probability=HK_R_probability';
    HK_R_vector(:,3)=HK_R_probability;
    
    total_rate(2,1)=sum(HK_R_probability);
    total_rate(2,2)=0;
    
    HK_R_probability=HK_R_probability/(min(HK_R_probability))+s;
    
    HK_R_probability=cumsum(log(HK_R_probability).^pf);    
    
    HK_R_vector(:,1)=HK_R_probability/HK_R_probability(end);
    IL_R_vector=[];
    
    % Recombination    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hopping
    
    HK_H_vector=zeros(hkl,3);
    HK_H_vector(:,2)=OV_list{1};

    
    HK_H_probability=Ks*factor_h_HK*exp(-Eh_h./(K*T_HK(OV_list{1}))).*exp(p_h*(KH(OV_list{1})+2)/3.*E_HK(OV_list{1})./(K*T_HK(OV_list{1})));    
    HK_H_probability=HK_H_probability'; 
    HK_H_vector(:,3)=HK_H_probability;
    
    total_rate(3,1)=sum(HK_H_probability);
    total_rate(3,2)=0;
    
    HK_H_probability=HK_H_probability/(min(HK_H_probability))+s;
    
    HK_H_probability=cumsum(log(HK_H_probability).^pf);  
    
    HK_H_vector(:,1)=HK_H_probability/HK_H_probability(end);
    IL_H_vector=[];    
       
else

    HK_R_vector=zeros(hkl,3);
    HK_R_vector(:,2)=OV_list{1};
    IL_R_vector=zeros(ill,3);
    IL_R_vector(:,2)=OV_list{2};
    
    HK_R_probability=factor_r_HK*exp((T_HK(OV_list{1})-Initial_T)/10+3).*exp(-Er_h./(K*T_HK(OV_list{1})));
    IL_R_probability=factor_r_IL*exp((T_IL(OV_list{2})-Initial_T)/30+3).*exp(-Er_s./(K*T_IL(OV_list{2})));
    
    HK_R_probability=HK_R_probability';
    IL_R_probability=IL_R_probability';
    HK_R_vector(:,3)=HK_R_probability;
    IL_R_vector(:,3)=IL_R_probability;        
    
    total_rate(2,1)=sum(HK_R_probability);
    total_rate(2,2)=sum(IL_R_probability);
    
    HK_R_probability=HK_R_probability/(min(HK_R_probability))+s;
    IL_R_probability=IL_R_probability/(min(IL_R_probability))+s;
    
    HK_R_probability=cumsum(log(HK_R_probability).^pf);
    IL_R_probability=cumsum(log(IL_R_probability).^pf);
    
    HK_R_vector(:,1)=HK_R_probability/HK_R_probability(end);
    IL_R_vector(:,1)=IL_R_probability/IL_R_probability(end);
    
    % Recombination    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hopping
    
    HK_H_vector=zeros(hkl,3);
    HK_H_vector(:,2)=OV_list{1};
    IL_H_vector=zeros(ill,3);
    IL_H_vector(:,2)=OV_list{2};
    
    HK_H_probability=Ks*factor_h_HK*exp(-Eh_h./(K*T_HK(OV_list{1}))).*exp(p_h*(KH(OV_list{1})+2)/3.*E_HK(OV_list{1})./(K*T_HK(OV_list{1})));
    IL_H_probability=Ks*factor_h_IL*exp(-Eh_s./(K*T_IL(OV_list{2}))).*exp(p_s*(KS(OV_list{2})+2)/3.*E_IL(OV_list{2})./(K*T_IL(OV_list{2})));
    
    HK_H_probability=HK_H_probability';
    IL_H_probability=IL_H_probability';
    HK_H_vector(:,3)=HK_H_probability;
    IL_H_vector(:,3)=IL_H_probability;      
        
    total_rate(3,1)=sum(HK_H_probability);
    total_rate(3,2)=sum(IL_H_probability);
    
    HK_H_probability=HK_H_probability/(min(HK_H_probability))+s;
    IL_H_probability=IL_H_probability/(min(IL_H_probability))+s; 
    
    HK_H_probability=cumsum(log(HK_H_probability).^pf);
    IL_H_probability=cumsum(log(IL_H_probability).^pf); 
    
    HK_H_vector(:,1)=HK_H_probability/HK_H_probability(end);
    IL_H_vector(:,1)=IL_H_probability/IL_H_probability(end);    
    
end

Probability(1,1)={HK_G_vector};
Probability(1,2)={IL_G_vector};
Probability(2,1)={HK_R_vector};
Probability(2,2)={IL_R_vector};
Probability(3,1)={HK_H_vector};
Probability(3,2)={IL_H_vector};

%% MCMC Estimation Parameters

if c1==1&&c2==1
    
    total_rate(:,3:4)=1;
    MC=sum(total_rate(:,1:2),1);
    MC=log(MC);
    MC=MC-min(MC)+s1;
    
    MC=log(MC)+s1;
    MC=log(MC)+s1;   
    MC=cumsum(log(MC));
    total_rate(:,5)=MC(1)/MC(2);
    
%   ind=min(find(abs(total_rate(:,4)-p)==min(abs(total_rate(:,4)-p))));

elseif c1==1&&c2==0
    
    KM_rate=total_rate(:,2);
    KM_rate=log(KM_rate);
    KM_rate=bsxfun(@minus,KM_rate,min(KM_rate,[],1))+1;
    KM_rate=cumsum(KM_rate,1);
    KM_rate=bsxfun(@rdivide,KM_rate,KM_rate(3,:));
    total_rate(:,4)=KM_rate;
    total_rate(:,3)=1;
    
    MC=sum(total_rate(:,1:2),1);
    MC=log(MC);
    MC=MC-min(MC)+s1;
    MC=log(MC)+s1;
    MC=log(MC)+s1;
    MC=cumsum(log(MC));
    total_rate(:,5)=MC(1)/MC(2);
    
elseif c2==1&&c1==0
    
    KM_rate=total_rate(:,1);
    KM_rate=log(KM_rate);
    KM_rate=bsxfun(@minus,KM_rate,min(KM_rate,[],1))+1;
    KM_rate=cumsum(KM_rate,1);
    KM_rate=bsxfun(@rdivide,KM_rate,KM_rate(3,:));
    total_rate(:,3)=KM_rate;
    total_rate(:,4)=1;  
    
    MC=sum(total_rate(:,1:2),1);
    MC=log(MC);
    MC=MC-min(MC)+s1;
    MC=log(MC)+s1;
    MC=log(MC)+s1;
    MC=cumsum(log(MC));
    total_rate(:,5)=MC(1)/MC(2);
    
else
    
    KM_rate=total_rate(:,1:2);
    KM_rate=log(KM_rate);
    KM_rate=bsxfun(@minus,KM_rate,min(KM_rate,[],1))+1;
    KM_rate=cumsum(KM_rate,1);
    KM_rate=bsxfun(@rdivide,KM_rate,KM_rate(3,:));
    total_rate(:,3:4)=KM_rate;
    
    MC=sum(total_rate(:,1:2),1);
    MC=log(MC);
    MC=MC-min(MC)+s1;
    MC=log(MC)+s1;
    MC=log(MC)+s1;
    MC=cumsum(log(MC));
    total_rate(:,5)=MC(1)/MC(2);
    
    tune=total_rate(:,1);
    tune=log(tune);
    tune=tune-min(tune)+4;
    tune=cumsum(tune);
    tune=tune/tune(end);
    total_rate(:,3)=tune;
    
end

%% TAT Model and C matrix Update
%C(:,:,1:T_hk)=C_hk*exp(-Ec_h./(K*T(:,:,1:T_hk)));
%toc
end