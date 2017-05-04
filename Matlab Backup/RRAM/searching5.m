function [podes,S_rate,sped,R_rate]=searching5(position,clas,poC,th_direction,Idensity,Incretemp)
global C Diffusion O n1 n2 K T e Z Qs Da omega Module_S N Unitcon Res tcr Initial_T pr th_ex Ks
%% Initial
t=max(4*n1+1,2*n2-1);
c=0:1:t;
d=Diffusion(clas(1)).dir(clas(2),:);
p=position;
C_temp=T(p(1),p(2),p(3))+Incretemp;
pC=poC;
h_flow=(-th_direction(1)*d(1)-th_direction(2)*d(2)*sqrt(2)-th_direction(3)*d(3)*sqrt(2))/(sqrt(d(1)^2+2*d(2)^2+2*d(3)^2)*20);
%% Find empty location for possible motion
if p(1)<=2||p(1)>=4*n1||p(2)<=2||p(2)>=4*n1
    podes=[0,1,1];%boundary conditions
    S_rate=0;sped=0;R_rate=zeros(1,4);
    return
else
    if d(1)>0
        a1=find(4*n1+1-(p(1)+d(1)*c)==2);
        a2=find(4*n1+1-(p(1)+d(1)*c)==3);
        a3=find(p(1)-d(1)*c==2);
        a4=find(p(1)-d(1)*c==3);
        ax=[a1 a2];ay=[a3 a4];
    elseif d(1)<0
        a1=find(p(1)+d(1)*c==2);
        a2=find(p(1)+d(1)*c==3);
        a3=find(4*n1+1-(p(1)-d(1)*c)==2);
        a4=find(4*n1+1-(p(1)-d(1)*c)==3);
        ax=[a1 a2];ay=[a3 a4];
    else
        ax=inf;ay=inf;
    end
    if d(2)>0
        b1=find(4*n1+1-(p(2)+d(2)*c)==2);
        b2=find(4*n1+1-(p(2)+d(2)*c)==3);
        b3=find(p(2)-d(2)*c==2);
        b4=find(p(2)-d(2)*c==3);
        bx=[b1 b2];by=[b3 b4];
    elseif d(2)<0
        b1=find(p(2)+d(2)*c==2);
        b2=find(p(2)+d(2)*c==3);
        b3=find(4*n1+1-(p(2)-d(2)*c)==2);
        b4=find(4*n1+1-(p(2)-d(2)*c)==3);
        bx=[b1 b2];by=[b3 b4];
    else
        bx=inf;by=inf;
    end
    if d(3)>0
        f1=find(2*n2-1-(p(3)+d(3)*c)==2);
        f2=find(2*n2-1-(p(3)+d(3)*c)==3);
        f3=find(p(3)-d(3)*c==2);
        f4=find(p(3)-d(3)*c==3);
        fx=[f1 f2];fy=[f3 f4];
    else
        fx=inf;fy=inf;
    end
    l=min([ax bx fx]);m=min([ay by fy]);
    if isempty(ay)||isempty(by)||isempty(fy)
        j=0:1:l-1;
        ind1=sub2ind(size(O),p(1)+d(1)*j,p(2)+d(2)*j,p(3)+d(3)*j);
        c1=O(ind1);c2=C(ind1);
        if ~isempty(find(c2==2, 1))
            fprintf('Occupation error');
            fprintf('%d\n',p,d,find(c2==2,1));
            return
        end
        k=find(c1==1,1);nu=length(find(c1==1));o=2;
    else
        j=0:1:l-1;i=0:1:m-1;
        ind1=sub2ind(size(O),p(1)+d(1)*j,p(2)+d(2)*j,p(3)+d(3)*j);
        ind2=sub2ind(size(O),p(1)-d(1)*i,p(2)-d(2)*i,p(3)-d(3)*i);
        c1=O(ind1);c2=C(ind1);c3=O(ind2);
        if ~isempty(find(c2==2, 1))
            fprintf('Occupation error');
            fprintf('%d\n',p,d,find(c2==2,1));
            return
        end
        k=find(c1==1,1);nu=length(find(c1==1));o=find(c3==1,1);
        if isempty(o)
            o=m;
        end
    end
%% k-2:number of occupied positions(forward), nu:number of non-occupied position(forward),o-2:number of occupied position(backward)
%% l-1:number of total position(forward), m-1:number of total position (backward).
%% Decide the rates
    if d(3)==0%% plane diffusion condition
        if  isempty(k)
            podes=[0,1,1];%boundary conditions
            S_rate=0;sped=0;R_rate=zeros(1,4);
            return
        else
%% Concentration at destination position
            podes=[p(1)+d(1)*(k-1),p(2)+d(2)*(k-1),p(3)+d(3)*(k-1)];
            if (podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1)&&(podes(3)~=1)
                sample=O(podes(1)-1:1:podes(1)+1,podes(2)-1:1:podes(2)-1,podes(3)-1:1:podes(3)+1);
                n=length(find(sample==3));
                podesC=(n/27)*(125/18);
            elseif (podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1)&&(podes(3)==1)
                sample=O(podes(1)-1:1:podes(1)+1,podes(2)-1:1:podes(2)+1,podes(3):1:podes(3)+2);
                n=length(find(sample==3));
                podesC=(n/27)*(125/18);
            elseif (podes(3)==2)&&(~(podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1))
                sample=O(podes(1)-2:1:podes(1)+2,podes(2)-2:1:podes(2)+2,podes(3)-1:1:podes(3)+3);
                n=length(find(sample==3));
                podesC=n/18;
            elseif (podes(3)==1)&&(~(podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1))
                sample=O(podes(1)-2:1:podes(1)+2,podes(2)-2:1:podes(2)+2,podes(3):1:podes(3)+4);
                n=length(find(sample==3));
                podesC=n/18;
            else
                sample=O(podes(1)-2:1:podes(1)+2,podes(2)-2:1:podes(2)+2,podes(3)-2:1:podes(3)+2);
                n=length(find(sample==3));
                podesC=n/18;
            end
            if ((2*n1+1-p(1))*d(1)+(2*n1+1-p(2))*d(2))>0 %Decide the path condition
                [~,v]=min(abs((d(2)/d(1))*(p(1)+d(1)*j-2*n1-1)+2*n1+1-p(2)-d(2)*j));
                k=v(1)+1;
                R_em=-double((pC*N/omega)*(o/((o+k-2)))*Res*(1+tcr*(C_temp-Initial_T))*(e*Z*Idensity/(K*C_temp))*Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp)));% Electromigration rate
                R_df=double(Da*exp(-(0.73*e)/(K*C_temp))*((pC*N/omega)-((podesC*N/omega)*(k-2)+(pC*N/omega))/(k-1))*(1/((k-1)*Unitcon*Diffusion(clas(1)).dis)));% Diffusion rate
                R_st=double((Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*(omega/N)*(th_ex*Module_S*h_flow/(1-2*pr))));% Stress diffusion rate
                R_th=double(Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*Qs*(h_flow/C_temp));% Thermal diffusion rate
                S_rate1=R_em+R_df+R_st+R_th;
                if S_rate1>0
                    S_rate=S_rate1;
                else
                    S_rate=0;
                    podes=[0,1,1];sped=0;R_rate=zeros(1,4);
                    return
                end
            elseif ((2*n1+1-p(1))*d(1)+(2*n1+1-p(2))*d(2))<0
                [~,v]=min(abs((d(2)/d(1))*(p(1)-d(1)*i-2*n1-1)+2*n1+1-p(2)+d(2)*i));
                o=v(1)+1;
                R_em=double((pC*N/omega)*(o/((o+k-2)))*Res*(1+tcr*(C_temp-Initial_T))*(e*Z*Idensity/(K*C_temp))*Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp)));% Electromigration rate
                R_df=double(Da*exp(-(0.73*e)/(K*C_temp))*((pC*N/omega)-((podesC*N/omega)*(k-2)+(pC*N/omega))/(k-1))*(1/((k-1)*Unitcon*Diffusion(clas(1)).dis)));% Diffusion rate
                R_st=double((Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*(omega/N)*(th_ex*Module_S*h_flow/(1-2*pr))));% Stress diffusion rate
                R_th=double(Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*Qs*(h_flow/C_temp));% Thermal diffusion rate
                S_rate1=R_em+R_df+R_st+R_th;
                if S_rate1>0
                    S_rate=S_rate1;
                else
                    S_rate=0;
                    podes=[0,1,1];sped=0;R_rate=zeros(1,4);
                    return
                end
            else
                o=2;
                R_em=double((pC*N/omega)*(o/((o+k-2)))*Res*(1+tcr*(C_temp-Initial_T))*(e*Z*Idensity/(K*C_temp))*Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp)));% Electromigration rate
                R_df=double(Da*exp(-(0.73*e)/(K*C_temp))*((pC*N/omega)-((podesC*N/omega)*(k-2)+(pC*N/omega))/(k-1))*(1/((k-1)*Unitcon*Diffusion(clas(1)).dis)));% Diffusion rate
                R_st=double((Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*(omega/N)*(th_ex*Module_S*h_flow/(1-2*pr))));% Stress diffusion rate
                R_th=double(Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*Qs*(h_flow/C_temp));% Thermal diffusion rate
                S_rate1=R_em+R_df+R_st+R_th;
                if S_rate1>0
                    S_rate=S_rate1;
                else
                    S_rate=0;
                    podes=[0,1,1];sped=0;R_rate=zeros(1,4);
                    return
                end
            end
        end
    else %% non-plane diffusion condition
        if  isempty(k)
            podes=[0,1,1];%boundary conditions
            S_rate=0;sped=0;R_rate=zeros(1,4);
        else
%% Concentration at destination position            
            podes=[p(1)+d(1)*(k-1),p(2)+d(2)*(k-1),p(3)+d(3)*(k-1)];
            if (podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1)&&(podes(3)~=1)
                sample=O(podes(1)-1:1:podes(1)+1,podes(2)-1:1:podes(2)-1,podes(3)-1:1:podes(3)+1);
                n=length(find(sample==3));
                podesC=n/27;
            elseif (podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1)&&(podes(3)==1)
                sample=O(podes(1)-1:1:podes(1)+1,podes(2)-1:1:podes(2)+1,podes(3):1:podes(3)+2);
                n=length(find(sample==3));
                podesC=n/27;
            elseif (podes(3)==2)&&(~(podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1))
                sample=O(podes(1)-2:1:podes(1)+2,podes(2)-2:1:podes(2)+2,podes(3)-1:1:podes(3)+3);
                n=length(find(sample==3));
                podesC=n/125;
            elseif (podes(3)==1)&&(~(podes(1)<3||podes(1)>4*n1-1||podes(2)<3||podes(2)>4*n1-1))
                sample=O(podes(1)-2:1:podes(1)+2,podes(2)-2:1:podes(2)+2,podes(3):1:podes(3)+4);
                n=length(find(sample==3));
                podesC=n/125;
            else
                sample=O(podes(1)-2:1:podes(1)+2,podes(2)-2:1:podes(2)+2,podes(3)-2:1:podes(3)+2);
                n=length(find(sample==3));
                podesC=n/125;
            end
            R_em=double((pC*N/omega)*(o/((o+k-2)))*Res*(1+tcr*(C_temp-Initial_T))*(e*Z*Idensity/(K*C_temp))*Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp)));% Electromigration rate
            R_df=double(Da*exp(-(0.73*e)/(K*C_temp))*((pC*N/omega)-((podesC*N/omega)*(o+k-2)+(pC*N/omega))/(o+k-1))*(1/((k-1)*Unitcon*Diffusion(clas(1)).dis)));% Diffusion rate
            R_st=double((Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*(omega/N)*(th_ex*Module_S*h_flow/(1-2*pr))));% Stress diffusion rate
            R_th=double(Da*exp(-(Diffusion(clas(1)).aEnergy*e)/(K*C_temp))*((pC*N/omega)/(K*C_temp))*Qs*(h_flow/C_temp));% Thermal diffusion rate
            S_rate1=R_em+R_df+R_st+R_th;
            if S_rate1>0
                S_rate=S_rate1;
            else
                S_rate=0;
                podes=[0,1,1];sped=0;R_rate=zeros(1,4);
                return
            end
        end
    end
end
R_rate=[R_em,R_df,R_st,R_th];
sped=1/((S_rate*pi*sqrt(2)*(0.2*n1)^2*Unitcon^2)*Ks);% Speed of migration
end

