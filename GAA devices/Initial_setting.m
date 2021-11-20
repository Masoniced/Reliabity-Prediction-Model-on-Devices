function [P,C_diag,IVR_lable,R_lable,Location_lable,location_conduct,c,Kth,Thermal_contribution,D_K,Pre_check_value,R]=Initial_setting(n1,n2,n3,T_hk,C_hk,C_il,Kth_HK,Kth_IL,PIT_p)
global kh_gb ks_p kh ks_v
%tic
IVR_lable=zeros(2*n1-1,2*n2-1,2*n3-1);
IVR_lable(1:2:2*n1-1,1:2:2*n2-1,1:2:2*n3-1)=1;
IVR_lable(1:2:2*n1-1,2:2:2*n2-2,1:2:2*n3-1)=-1;
IVR_lable(2:2:2*n1-2,1:2:2*n2-1,1:2:2*n3-1)=-1;
IVR_lable(1:2:2*n1-1,1:2:2*n2-1,2:2:2*n3-2)=-1;


line_number=3*n1*n2*n3-n1*n2-n2*n3-n1*n3-(2*n1*n2-n1-n2)*2;
node_number=n1*n2*n3-n1*n2*2+2;
count=0;
D_K=zeros(n1,n2,n3-2);


Thermal_contribution=ones(line_number+1,1);
r_vec=find(IVR_lable==-1);
R_lable=zeros(line_number,12);
c=zeros(line_number,1);
Kth=zeros(line_number,1);
p1y=zeros(line_number,1);
p2y=zeros(line_number,1);
p1z=zeros(line_number,1);
p2z=zeros(line_number,1);
for i=1:length(r_vec)
    [id1,id2,id3]=ind2sub([2*n1-1,2*n2-1,2*n3-1],r_vec(i));
    if (id3==1)||(id3==2*n3-1)
        continue
    elseif mod(id1,2)==0
%         if (IVR_lable(id1-1,id2,id3)~=1)||(IVR_lable(id1+1,id2,id3)~=1)
%             fprintf('intial lable error_id1')
%             fprintf('%d\n',id1,id2,id3)
%             break
%         end
        count=count+1;
        if id3<=2*T_hk+2
            c(count)=C_hk;
            Kth(count)=Kth_HK;
            R_lable(count,:)=[count,r_vec(i),id1,id2,id3,id1/2,(id2+1)/2,(id3-1)/2,id1/2+1,(id2+1)/2,(id3-1)/2,1];
        else
            c(count)=C_il;
            Kth(count)=Kth_IL;
            R_lable(count,:)=[count,r_vec(i),id1,id2,id3,id1/2,(id2+1)/2,(id3-1)/2,id1/2+1,(id2+1)/2,(id3-1)/2,2];
        end
        index1=sub2ind([n1,n2,n3-2],id1/2,(id2+1)/2,(id3+1)/2-1);
        index2=sub2ind([n1,n2,n3-2],(id1+2)/2,(id2+1)/2,(id3+1)/2-1);
        p1y(count)=index1+2;p1z(count)=-1;
        p2y(count)=index2+2;p2z(count)=1;
    elseif mod(id2,2)==0
%         if (IVR_lable(id1,id2-1,id3)~=1)||(IVR_lable(id1,id2+1,id3)~=1)
%             fprintf('intial lable error_id2')
%             fprintf('%d\n',id1,id2,id3)
%             break
%         end
        count=count+1;
        if id3<=2*T_hk+2
            c(count)=C_hk;
            Kth(count)=Kth_HK;            
            R_lable(count,:)=[count,r_vec(i),id1,id2,id3,(id1+1)/2,id2/2,(id3-1)/2,(id1+1)/2,id2/2+1,(id3-1)/2,1];
        else
            c(count)=C_il;
            Kth(count)=Kth_IL;
            R_lable(count,:)=[count,r_vec(i),id1,id2,id3,(id1+1)/2,id2/2,(id3-1)/2,(id1+1)/2,id2/2+1,(id3-1)/2,2];
        end        
        index1=sub2ind([n1,n2,n3-2],(id1+1)/2,id2/2,(id3+1)/2-1);
        index2=sub2ind([n1,n2,n3-2],(id1+1)/2,(id2+2)/2,(id3+1)/2-1);
        p1y(count)=index1+2;p1z(count)=-1;
        p2y(count)=index2+2;p2z(count)=1;
    elseif mod(id3,2)==0&&id3~=2&&id3~=2*n3-2
%         if (IVR_lable(id1,id2,id3-1)~=1)||(IVR_lable(id1,id2,id3+1)~=1)
%             fprintf('intial lable error_id3')
%             fprintf('%d\n',id1,id2,id3)
%             break
%         end
        count=count+1;
        if id3<=2*T_hk+2
            c(count)=C_hk;
            Kth(count)=Kth_HK;
            R_lable(count,:)=[count,r_vec(i),id1,id2,id3,(id1+1)/2,(id2+1)/2,id3/2-1,(id1+1)/2,(id2+1)/2,id3/2,1];
        else
            c(count)=C_il;
            Kth(count)=Kth_IL;
            R_lable(count,:)=[count,r_vec(i),id1,id2,id3,(id1+1)/2,(id2+1)/2,id3/2-1,(id1+1)/2,(id2+1)/2,id3/2,2];
        end         
        index1=sub2ind([n1,n2,n3-2],(id1+1)/2,(id2+1)/2,id3/2-1);
        index2=sub2ind([n1,n2,n3-2],(id1+1)/2,(id2+1)/2,(id3+2)/2-1);        
        p1y(count)=index1+2;p1z(count)=-1;
        p2y(count)=index2+2;p2z(count)=1;
    elseif id3==2
%         if (IVR_lable(id1,id2,id3-1)~=1)||(IVR_lable(id1,id2,id3+1)~=1)
%             fprintf('intial lable error_id3')
%             fprintf('%d\n',id1,id2,id3)
%             break
%         end        
        count=count+1;
        R_lable(count,:)=[count,r_vec(i),id1,id2,id3,(id1+1)/2,(id2+1)/2,id3/2-1,(id1+1)/2,(id2+1)/2,id3/2,1];
        c(count)=C_hk;
        Kth(count)=Kth_HK;
        index1=-1;
        index2=sub2ind([n1,n2,n3-2],(id1+1)/2,(id2+1)/2,(id3)/2);     
        p1y(count)=index1+2;p1z(count)=-1;
        p2y(count)=index2+2;p2z(count)=1;
    elseif id3==2*n3-2
%         if (IVR_lable(id1,id2,id3-1)~=1)||(IVR_lable(id1,id2,id3+1)~=1)
%             fprintf('intial lable error_id3')
%             fprintf('%d\n',id1,id2,id3)
%             break
%         end        
        count=count+1;
        R_lable(count,:)=[count,r_vec(i),id1,id2,id3,(id1+1)/2,(id2+1)/2,id3/2-1,(id1+1)/2,(id2+1)/2,id3/2,2];
        c(count)=C_il;
        Kth(count)=Kth_IL;
        index1=sub2ind([n1,n2,n3-2],(id1+1)/2,(id2+1)/2,id3/2-1);
        index2=0;      
        p1y(count)=index1+2;p1z(count)=-1;
        p2y(count)=index2+2;p2z(count)=1;
    end
end


C_diag=sparse(1:line_number,1:line_number,c);
P=sparse([1:line_number,1:line_number],[p1y,p2y],[p1z,p2z]);

Location_lable=zeros(n1*n2*(n3-2),6);
location_conduct=cell(n1*n2*(n3-2),1);

[R,Pre_check_value] = RPG_r(0,PIT_p,n1,n2,n3-2,T_hk);

for i=1:n1*n2*(n3-2)
    [ax,ay,az]=ind2sub([n1,n2,n3-2],i);
    p=[2*ax-1,2*ay-1,2*az+1];
    a1=2*n1-1;
    a2=2*n2-1;
    a3=2*n3-1;
    if R(i)==1
        if i<=n1*n2*T_hk
            D_K(i)=kh_gb;
        else
            D_K(i)=ks_p;
        end
    else
         if i<=n1*n2*T_hk
            D_K(i)=kh;
        else
            D_K(i)=ks_v;
         end
    end
    if p(1)==1&&p(2)==1&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2);  % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+(p(2))*a1+(p(3)-1)*a1*a2);% p(1)+(p(2))*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,line_number+1,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2)];

    elseif p(1)==1&&p(2)==a2&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2); % p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,line_number+1,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2)];

    elseif p(1)==a1&&p(2)==1&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+(p(2))*a1+(p(3)-1)*a1*a2); % p(1)+(p(2))*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,line_number+1,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2)];

    elseif p(1)==a1&&p(2)==a2&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2); % p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,line_number+1,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2)];

    %% Sides
    elseif p(1)==1&&p(2)~=1&&p(2)~=a2&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2); % p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2
        element5=find(R_lable(:,2)==p(1)+(p(2))*a1+(p(3)-1)*a1*a2); % p(1)+(p(2))*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,element5,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2),R_lable(element5,2)];

    elseif p(1)==a1&&p(2)~=1&&p(2)~=a2&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2); % p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2
        element5=find(R_lable(:,2)==p(1)+(p(2))*a1+(p(3)-1)*a1*a2); % p(1)+(p(2))*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,element5,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2),R_lable(element5,2)];

    elseif p(2)==1&&p(1)~=1&&p(1)~=a1&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element5=find(R_lable(:,2)==p(1)+(p(2))*a1+(p(3)-1)*a1*a2); % p(1)+(p(2))*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,element5,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2),R_lable(element5,2)];

    elseif p(2)==a2&&p(1)~=1&&p(1)~=a1&&p(3)~=1&&p(3)~=a3

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element5=find(R_lable(:,2)==p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2); % p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,element5,line_number+1];
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2),R_lable(element5,2)];

    %% Normal
    else

        element1=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2); % p(1)+(p(2)-1)*a1+(p(3)-2)*a1*a2
        element2=find(R_lable(:,2)==p(1)+(p(2)-1)*a1+(p(3))*a1*a2); % p(1)+(p(2)-1)*a1+(p(3))*a1*a2
        element3=find(R_lable(:,2)==p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)-1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element4=find(R_lable(:,2)==p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2); % p(1)+1+(p(2)-1)*a1+(p(3)-1)*a1*a2
        element5=find(R_lable(:,2)==p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2); % p(1)+(p(2)-2)*a1+(p(3)-1)*a1*a2
        element6=find(R_lable(:,2)==p(1)+(p(2))*a1+(p(3)-1)*a1*a2); % p(1)+(p(2))*a1+(p(3)-1)*a1*a2
        list=[element1,element2,element3,element4,element5,element6]; 
        list_c=[R_lable(element1,2),R_lable(element2,2),R_lable(element3,2),R_lable(element4,2),R_lable(element5,2),R_lable(element6,2)];

    end
    Location_lable(i,:)=list;
    location_conduct(i)={list_c};
    
end

%toc
end
