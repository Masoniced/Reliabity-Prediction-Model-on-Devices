%% Defination of Sampling region 
function [k]=sample(p,a1,a2,a3,value) % sub-function
global M 

%% For Add+Remove Clustering
kickoff=[p(1)-1+(p(2)-2)*a1+(p(3)-2)*a1*a2,p(1)-1+(p(2))*a1+(p(3)-2)*a1*a2,p(1)+1+(p(2)-2)*a1+(p(3)-2)*a1*a2,p(1)+1+(p(2))*a1+(p(3)-2)*a1*a2];
kickoff=[kickoff,p(1)-1+(p(2)-2)*a1+(p(3))*a1*a2,p(1)-1+(p(2))*a1+(p(3))*a1*a2,p(1)+1+(p(2)-2)*a1+(p(3))*a1*a2,p(1)+1+(p(2))*a1+(p(3))*a1*a2];
kickoff=sort(kickoff);
if value==1
%% Coners    

    if p(1)==1&&p(2)==1&&p(3)==1
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2):p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2):p(2)+1,p(3):p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==1&&p(2)==1&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2):p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2):p(2)+1,p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==1&&p(2)==a2&&p(3)==1
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2),p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2)-1:p(2),p(3):p(3)+1)==1);
        k(k==number)=[]; 
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==1&&p(2)==a2&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2),p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2)-1:p(2),p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(2)==1&&p(3)==1
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2):p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2):p(2)+1,p(3):p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(2)==a2&&p(3)==1
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2),p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2)-1:p(2),p(3):p(3)+1)==1);
        k(k==number)=[]; 
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(2)==1&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2):p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2):p(2)+1,p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(2)==a2&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2),p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2)-1:p(2),p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    %% Lines
    elseif p(1)==1&&p(2)==1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2):p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2):p(2)+1,p(3)-1:p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==1&&p(2)==a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2),p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2)-1:p(2),p(3)-1:p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(2)==1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2):p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2):p(2)+1,p(3)-1:p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(2)==a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2),p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2)-1:p(2),p(3)-1:p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==1&&p(3)==1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2)-1:p(2)+1,p(3):p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==1&&p(3)==a3&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(3)==1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2)-1:p(2)+1,p(3):p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(3)==a3&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2)-1:p(2)+1,p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(2)==1&&p(3)==1&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2):p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2):p(2)+1,p(3):p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(2)==1&&p(3)==a3&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2):p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2):p(2)+1,p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(2)==a2&&p(3)==1&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2),p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2)-1:p(2),p(3):p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(2)==a2&&p(3)==a3&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2),p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2)-1:p(2),p(3)-1:p(3))==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    %% Sides
    elseif p(1)==1&&p(2)~=1&&p(2)~=a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1):p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3)+1)==1);
        k(k==number)=[]; 
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(1)==a1&&p(2)~=1&&p(2)~=a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1),p(2)-1:p(2)+1,p(3)-1:p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(2)==1&&p(1)~=1&&p(1)~=a1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2):p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2):p(2)+1,p(3)-1:p(3)+1)==1);
        k(k==number)=[]; 
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(2)==a2&&p(1)~=1&&p(1)~=a1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2),p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2)-1:p(2),p(3)-1:p(3)+1)==1);
        k(k==number)=[]; 
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(3)==1&&p(1)~=1&&p(1)~=a1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3):p(3)+1)==1);
        k(k==number)=[]; 
        ind=ismembc(k,kickoff);
        k(ind)=[];
    elseif p(3)==a3&&p(1)~=1&&p(1)~=a1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3))==1);
        k(k==number)=[]; 
        ind=ismembc(k,kickoff);
        k(ind)=[];
    %% Normal
    else
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2(M(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3)+1)==1);
        k(k==number)=[];
        ind=ismembc(k,kickoff);
        k(ind)=[];
    end
    
else
%% For Deciede_hopping    

    if p(1)==1&&p(2)==1&&p(3)==1
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2):p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==1&&p(2)==1&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2):p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==1&&p(2)==a2&&p(3)==1
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2),p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];    
    elseif p(1)==1&&p(2)==a2&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2),p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==a1&&p(2)==1&&p(3)==1
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2):p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==a1&&p(2)==a2&&p(3)==1
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2),p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[]; 
    elseif p(1)==a1&&p(2)==1&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2):p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==a1&&p(2)==a2&&p(3)==a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2),p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    %% Lines
    elseif p(1)==1&&p(2)==1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2):p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==1&&p(2)==a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2),p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==a1&&p(2)==1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2):p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==a1&&p(2)==a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2),p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==1&&p(3)==1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==1&&p(3)==a3&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==a1&&p(3)==1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(1)==a1&&p(3)==a3&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(2)==1&&p(3)==1&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2):p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(2)==1&&p(3)==a3&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2):p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(2)==a2&&p(3)==1&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2),p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(2)==a2&&p(3)==a3&&p(1)~=1&&p(1)~=a1
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2),p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    %% Sides
    elseif p(1)==1&&p(2)~=1&&p(2)~=a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1):p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];    
    elseif p(1)==a1&&p(2)~=1&&p(2)~=a2&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1),p(2)-1:p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    elseif p(2)==1&&p(1)~=1&&p(1)~=a1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2):p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[]; 
    elseif p(2)==a2&&p(1)~=1&&p(1)~=a1&&p(3)~=1&&p(3)~=a3
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2),p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[]; 
    elseif p(3)==1&&p(1)~=1&&p(1)~=a1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3):p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[]; 
    elseif p(3)==a3&&p(1)~=1&&p(1)~=a1&&p(2)~=1&&p(2)~=a2
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3));
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[]; 
    %% Normal
    else
        [X,Y,Z]=meshgrid(p(1)-1:p(1)+1,p(2)-1:p(2)+1,p(3)-1:p(3)+1);
        s2=X+(Y-1)*a1+(Z-1)*a1*a2;
        s2=permute(s2,[2,1,3]);
        number=p(1)+(p(2)-1)*a1+(p(3)-1)*a1*a2;
        k=s2;
        k(k==number)=[];
    end    
end
end