%% TAT Model and C matrix Update
function [C_diag,s]=update_C(V,T_v,VI,scale)
global  n1 n2 n3 R_lable Ks K aphla_C Conduct E_con alph tunneling_lable e M
%tic
% R_label (line_number*12)
% [1: Resistor number; 2: Resistor Position Index; 3: Resistor Position X;
% 4: Resistor Position Y; 5: Resistor Position Z; 6: Connection point1
% Position X; 7: Connection point1 Position Y; 8: Connection point1
% Position Z; 9: Connection point2 Position X; 10: Connection point2
% Position Y; 11: Connection point2 Position Z; 12: Label for HK or IL;
if nargin==4
    num=log(scale)/log(log(400+3));
else
    num=1.6; % original: 1.2
end
factor=10;
damping=1e-8;
s_const=log(sum(M(:))+3)^num;
% s_const=1;
leng=(3*n1*n2*n3-n1*n2-n2*n3-n1*n3-(2*n1*n2-n1-n2)*2);
s=zeros(leng,1);
for i=1:leng
    if R_lable(i,5)==2
        s(i)=Conduct(R_lable(i,12))*exp(-E_con(R_lable(i,12))/(K*T_v(i)))*s_const....
            +damping*tunneling_lable(R_lable(i,3),R_lable(i,4),R_lable(i,5))*(e^2/(aphla_C*K*T_v(i)))....
            *Ks*exp(-2*aphla_C/alph)*exp(-e*abs((V(R_lable(i,9),R_lable(i,10),R_lable(i,11))-VI)/factor)/(K*T_v(i)));
    elseif R_lable(i,5)==2*n3-2
        s(i)=Conduct(R_lable(i,12))*exp(-E_con(R_lable(i,12))/(K*T_v(i)))*s_const....
            +damping*tunneling_lable(R_lable(i,3),R_lable(i,4),R_lable(i,5))*(e^2/(aphla_C*K*T_v(i)))....
            *Ks*exp(-2*aphla_C/alph)*exp(-e*(abs(V(R_lable(i,6),R_lable(i,7),R_lable(i,8))/factor))/(K*T_v(i)));
    else
        s(i)=Conduct(R_lable(i,12))*exp(-E_con(R_lable(i,12))/(K*T_v(i)))*s_const....
            +damping*tunneling_lable(R_lable(i,3),R_lable(i,4),R_lable(i,5))*(e^2/(aphla_C*K*T_v(i)))....
            *Ks*exp(-2*aphla_C/alph)*exp(-e*(abs(V(R_lable(i,6),R_lable(i,7),R_lable(i,8))-V(R_lable(i,9),R_lable(i,10),R_lable(i,11)))/factor)/(K*T_v(i)));
    end
end
C_diag=sparse(1:leng,1:leng,s);
%toc
end