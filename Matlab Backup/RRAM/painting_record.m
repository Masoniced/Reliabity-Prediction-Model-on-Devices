function [angle1,angle2,Z2,X2]=painting_record(O)
%% Painting
%% PLANE VIEW
uni_length=1.3575;% angstrom 1e-10m
Z=zeros(401,401);
for i=1:401
    L=O(:,i,:);
    L=reshape(L,401,401);
    Z=L+Z;
end
Z1=zeros(401,floor(401*uni_length));
Z1(:,1)=Z(:,1);
for i=2:floor(401*uni_length)
    x1=floor(i/uni_length);
    Z1(:,i)=Z(:,x1)+(i/uni_length-x1)*(Z(:,x1+1)-Z(:,x1));
end
Z2=zeros(floor(401*sqrt(2)*uni_length),floor(401*uni_length));
Z2(1,:)=Z1(1,:);
for i=2:floor(401*sqrt(2)*uni_length)
    x2=floor(i/(sqrt(2)*uni_length));
    Z2(i,:)=Z1(x2,:)+(i/(sqrt(2)*uni_length)-x2)*(Z1(x2+1,:)-Z1(x2,:));
end
Z2=Z2-401;
Z2=Z2';
Z2=Z2*(5.43)/4*sqrt(2);
kin=1;
for j=1:401;
    para=find(Z2(:,j),1,'last');
    if para>kin
        kin=para;
    end
end
bin1=(find(Z2(1,:),1,'first')+find(Z2(2,:),1,'first')+find(Z2(3,:),1,'first'))/3;
bin2=(find(Z2(1,:),1,'last')+find(Z2(2,:),1,'last')+find(Z2(2,:),1,'last'))/3;
tin1=(find(Z2(kin-3,:),1,'first')+find(Z2(kin-4,:),1,'first')+find(Z2(kin-5,:),1,'first'))/3;
tin2=(find(Z2(kin-3,:),1,'last')+find(Z2(kin-4,:),1,'last')+find(Z2(kin-5,:),1,'last'))/3;
angle1=180*atan((kin-7)/(tin1-bin1))/pi;angle2=180*atan((kin-7)/(bin2-tin2))/pi;
%% TOP VIEW
X=zeros(401,401);
for i=1:401
    for j=1:401
        if isempty(find(O(i,j,:)==3,1,'last'))
            X(i,j)=0;
        else
            X(i,j)=find(O(i,j,:)==3,1,'last');
        end
    end
end
X1=zeros(401,floor(401*sqrt(2)*uni_length));
for i=2:floor(401*sqrt(2)*uni_length)
    X1(:,1)=X(:,1);
    x3=floor(i/(sqrt(2)*uni_length));
    X1(:,i)=X(:,x3)+(i/(sqrt(2)*uni_length)-x3)*(X(:,x3+1)-X(:,x3));
end
X2=zeros(floor(401*sqrt(2)*uni_length),floor(401*sqrt(2)*uni_length));
for i=2:floor(401*sqrt(2)*uni_length)
    X2(1,:)=X1(1,:);
    x4=floor(i/(sqrt(2)*uni_length));
    X2(i,:)=X1(x4,:)+(i/(sqrt(2)*uni_length)-x4)*(X1(x4+1,:)-X1(x4,:));
end
X2=X2*(5.43)/4*sqrt(2);
end

