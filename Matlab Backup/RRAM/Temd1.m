function [dtemper]=Temd1(a,b,c,C_time,C_temp)
global n1 T Cur Unitcon Res MD_ns Cp_ns Initial_T tcr tk_ns 
Idensity=(Cur/(2*pi*Unitcon^2))*(1/(sqrt((2*(a-2*n1)^2+2*(b-2*n1)^2+(c)^2))+sqrt(1/2)*(0.2*n1+1))^2);% changed by position
if (T(a,b,c)~=300)
    dtemper=T(a,b,c);
else
    sod_t_time=12*Idensity^2*Res*(1+tcr*(C_temp-Initial_T))*(1/((sqrt((2*(a-2*n1)^2+2*(b-2*n1)^2+(c)^2))+sqrt(1/2)*(0.2*n1+1))^2*Unitcon^2))/((MD_ns)*(Cp_ns));
    dtemper=Initial_T+(Idensity^2*Res*(1+tcr*(C_temp-Initial_T))*(C_time)+sod_t_time*tk_ns*(C_time)^(2.5)/2)/((MD_ns)*(Cp_ns));
end
end

    