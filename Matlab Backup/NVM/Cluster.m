%% Dynamic Cluster Checking (Hoshen Kopelman A)
function [criteria,fb,decide_label_1]=Cluster(p) % main function
global Lable_mat M label Label_record
[a1,a2,a3]=size(M);
x=p(1);
y=p(2);
z=p(3);
u=sample(p,a1,a2,a3);
t1=length(u);
check_volum=[]; % exclude non-defined neighbours
for m=1:t1
    if Lable_mat(u(m))==0
            u(m)=inf;
    end
end
u(u==inf)=[];
o=length(u);
if t1-o>1
    fprintf('Sample error')
    fprintf('%d\n',p,u,o,t1)
    return
end
if p(3)==a3-1||p(3)==a3
    ci=1;
else
    ci=0;
end
if o==0 
    %% 0 neighbour case
    decide_label=label+1;
    Lable_mat(x,y,z)=decide_label;
    Label_record(decide_label,1)={[x,y,z]};
    Label_record(decide_label,2)={decide_label};
    label=decide_label;
elseif o==1 
    %% 1 neighbour case
    decide_label=Lable_mat(u);
    Lable_mat(x,y,z)=decide_label;
    Label_record(decide_label,1)={[Label_record{decide_label,1};x,y,z]};
    Label_record(decide_label,2)={decide_label};
else
    %% >1 neighbour case
    for j=1:o
        label_box=zeros(1,o);
        label_box(j)=Lable_mat(u(j));
    end
    [~,p2]=sort(label_box);
    decide_label=Lable_mat(u(p2(1)));
    Lable_mat(x,y,z)=decide_label;
    Label_record(decide_label,1)={[Label_record{decide_label,1};x,y,z]};
    for l=1:o
        Label_record(Lable_mat(u(l)),2)={decide_label};
    end
end
decide_label_1=decide_label;
%% Label combination
vector=[Label_record{:,2}];
check_vector=find(vector==decide_label);
w=length(check_vector);
for i=1:w
    check_volum=[check_volum;Label_record{check_vector(i),1}];
end
if ci==0
    fb=0;
else
    if (~isempty(find(check_volum(:,3)==a3-1, 1)))&&(~isempty(find(check_volum(:,3)==a3, 1)))
        fb=1;
    else
        fb=0;
    end
end   
if (~isempty(find(check_volum(:,3)==1, 1)))&&(~isempty(find(check_volum(:,3)==a3, 1)))
    criteria=1;
else
    criteria=0;
end  
end


