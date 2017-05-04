function [S_T]=smothing(T,x,y,z)
h=length(x);
T_sample=T;
S_T=zeros(401,401,401);
sig=3;
normal_kernel_3=zeros(3,3,3);
normal_kernel_5=zeros(5,5,5);
% normal kernel
for i=1:3
    for j=1:3
        for l=1:3
            normal_kernel_3(i,j,l)=(1/(2*pi*sig^2))*exp(-(abs(i-2)^2+abs(j-2)^2+abs(l-2)^2)/(2*sig^2));
        end
    end
end
for i=1:5
    for j=1:5
        for l=1:5
            normal_kernel_5(i,j,l)=(1/(2*pi*sig^2))*exp(-(abs(i-3)^2+abs(j-3)^2+abs(l-3)^2)/(2*sig^2));
        end
    end
end
samp_3=zeros(3,3,3);
samp_5=zeros(5,5,5);
for i=1:h
    if z(i)==1
        samp_3(:,:,1)=zeros(3,3);
        samp_3(:,:,2:1:3)=T_sample(x(i)-1:1:x(i)+1,y(i)-1:1:y(i)+1,z(i):1:z(i)+1);
        index=find(samp_3);
        t=length(index);
        weight=0;
        for j=1:t
            weight=weight+normal_kernel_3(index(t));
        end
        if weight==0
            fprintf('weight_error')
            fprintf('%d\n',x(i),y(i),z(i))
            return
        end
        S_T(x(i),y(i),z(i))=sum(sum(sum((samp_3.*normal_kernel_3))))/weight;
    elseif x(i)==1
        samp_3(1,:,:)=zeros(3,3);
        samp_3(2:1:3,:,:)=T_sample(x(i):1:x(i)+1,y(i)-1:1:y(i)+1,z(i)-1:1:z(i)+1);
        index=find(samp_3);
        t=length(index);
        weight=0;
        for j=1:t
            weight=weight+normal_kernel_3(index(t));
        end
        if weight==0
            fprintf('weight_error')
            fprintf('%d\n',x(i),y(i),z(i))
            return
        end
        S_T(x(i),y(i),z(i))=sum(sum(sum((samp_3.*normal_kernel_3))))/weight;
    elseif y(i)==1
        samp_3(:,1,:)=zeros(3,3);
        samp_3(:,2:1:3,:)=T_sample(x(i)-1:1:x(i)+1,y(i):1:y(i)+1,z(i)-1:1:z(i)+1);
        index=find(samp_3);
        t=length(index);
        weight=0;
        for j=1:t
            weight=weight+normal_kernel_3(index(t));
        end
        if weight==0
            fprintf('weight_error')
            fprintf('%d\n',x(i),y(i),z(i))
            return
        end
        S_T(x(i),y(i),z(i))=sum(sum(sum((samp_3.*normal_kernel_3))))/weight;
    elseif y(i)==2||x(i)==2
        samp_3=T_sample(x(i)-1:1:x(i)+1,y(i)-1:1:y(i)+1,z(i)-1:1:z(i)+1);
        index=find(samp_3);
        t=length(index);
        weight=0;
        for j=1:t
            weight=weight+normal_kernel_3(index(t));
        end
        if weight==0
            fprintf('weight_error')
            fprintf('%d\n',x(i),y(i),z(i))
            return
        end
        S_T(x(i),y(i),z(i))=sum(sum(sum((samp_3.*normal_kernel_3))))/weight;        
    elseif z(i)==2
        samp_5(:,:,1)=zeros(5,5);
        samp_5(:,:,2:1:5)=T_sample(x(i)-2:1:x(i)+2,y(i)-2:1:y(i)+2,z(i)-1:1:z(i)+2);
        index=find(samp_5);
        t=length(index);
        weight=0;
        for j=1:t
            weight=weight+normal_kernel_5(index(t));
        end
        if weight==0
            fprintf('weight_error')
            fprintf('%d\n',x(i),y(i),z(i))
            return
        end
        S_T(x(i),y(i),z(i))=sum(sum(sum((samp_5.*normal_kernel_5))))/weight;
    else
        samp_5=T_sample(x(i)-2:1:x(i)+2,y(i)-2:1:y(i)+2,z(i)-2:1:z(i)+2);
        index=find(samp_5);
        t=length(index);
        weight=0;
        for j=1:t
            weight=weight+normal_kernel_5(index(t));
        end
        if weight==0
            fprintf('weight_error')
            fprintf('%d\n',x(i),y(i),z(i))
            return
        end
        S_T(x(i),y(i),z(i))=sum(sum(sum((samp_5.*normal_kernel_5))))/weight;
    end
end
end
            
