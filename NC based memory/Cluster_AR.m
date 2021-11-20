%% Dynamic Cluster Checking (Hoshen Kopelman A)
function [criteria,fb,decide_label_1,check_volum_find]=Cluster_AR(p,level,n1,n2,n3) % main function
global Lable_mat label Label_record M

a1=n1;
a2=n2;
a3=n3-2;
x=p(1);
y=p(2);
z=p(3);

T_hk=4;

switch level
    case 0
        %% Add Component    
        u=sample(p,a1,a2,a3,1);
        t1=length(u);
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
        if p(3)>T_hk
            ci=1;
        else
            ci=0;
        end
        if o==0 
            %% 0 neighbour case
            decide_label=label+1;
            Lable_mat(x,y,z)=decide_label;
            Label_record(decide_label,1)={[x;y;z]};
            Label_record(decide_label,2)={x+(y-1)*a1+(z-1)*a1*a2};
            Label_record(decide_label,3)={decide_label};
            label=decide_label;
        elseif o==1 
            %% 1 neighbour case
            decide_label=Lable_mat(u);
            Lable_mat(x,y,z)=decide_label;
            Label_record(decide_label,1)={[[Label_record{decide_label,1}],[x;y;z]]};
            Label_record(decide_label,2)={[Label_record{decide_label,2},x+(y-1)*a1+(z-1)*a1*a2]};
            Label_record(decide_label,3)={decide_label};
        else
            %% >1 neighbour case
            u1=min(u);
            decide_label=Lable_mat(u1);
            u(u==u1)=[];
            Lable_mat(x,y,z)=decide_label;
            Label_record(decide_label,1)={[[Label_record{decide_label,1}],[x;y;z]]};
            Label_record(decide_label,2)={[Label_record{decide_label,2},x+(y-1)*a1+(z-1)*a1*a2]};
            box=sort([Label_record{Lable_mat(u),3}]);
            box=box([true,diff(box)>0]); % speed up the set-diff
            Label_record(ismembc([Label_record{:,3}],box),3)={decide_label};
        end
        decide_label_1=decide_label;
        %% Label combination
        check_volum=[Label_record{[Label_record{:,3}]==decide_label,1}];
        if ci==0
            fb=0;
        else
            if (~isempty(find(check_volum(3,:)==a3-T_hk+1, 1)))&&(~isempty(find(check_volum(3,:)==a3, 1)))
                fb=1;
            else
                fb=0;
            end
        end   
        if (~isempty(find(check_volum(3,:)==1, 1)))&&(~isempty(find(check_volum(3,:)==a3, 1)))
            criteria=1;
            check_volum_find=check_volum;
        else
            criteria=0;
            check_volum_find=[];
        end  
    case 1
        %% Remove Component
        cluster_label=Label_record{Lable_mat(x,y,z),3};
        Lable_mat(x,y,z)=0;
        M(x,y,z)=0;
        s_list=logical([Label_record{:,3}]==cluster_label);
        cluster_volum=[Label_record{s_list,1}]';
        cluster_volum2=[Label_record{s_list,2}]';
        inds=find(cluster_volum2==(x+(y-1)*a1+(z-1)*a1*a2));
        cluster_volum(inds,:)=[];
        cluster_volum2(inds)=[];
        Label_record(s_list,1)={[]};
        Label_record(s_list,2)={[]};
        decide_l=0;
        
        label_r={};
        if isempty(cluster_volum)
            num=0;
        else
            num=length(cluster_volum(:,1));
            if num>1    
                list=zeros(num,1);
                count=0;
                count1=0;
                criteria1=1;
                ss=1:num;
                while criteria1==1
                    vector=abs(bsxfun(@minus,cluster_volum,cluster_volum(count+1,:)));
                    value=vector(:,1).^2+vector(:,2).^2+vector(:,3).^2;
                    count=count+1;
                    ind=ss(value<=3);
                    if isempty(label_r)
                        u_o=[];
                    else
                        u_o=ismembc(ind,sort([label_r{:,1}]));% check if ind is already labeled in label_r
                    end
                    decide_l=decide_l+1;
                    con1=sum(u_o);
                    con2=length(ind);
                    if con2>con1&&con1~=0
                        count1=count1+1;
                        u_i=~u_o;
                        intersect=ind(u_o);
                        label_r(count1,1)={ind(u_i)};
                        label_r(count1,2)={decide_l};
                        list(ind(u_i))=count1;
                        r_box=sort([label_r{list(intersect),2}]);
                        r_box=r_box([true,diff(r_box)>0]);
                        label_r(ismembc([label_r{:,2}],r_box),2)={decide_l};
                    elseif con2==con1
                        r_box=sort([label_r{list(ind),2}]);
                        r_box=r_box([true,diff(r_box)>0]);
                        label_r(ismembc([label_r{:,2}],r_box),2)={decide_l};   
                    elseif con1==0
                        count1=count1+1;
                        label_r(count1,1)={ind};
                        label_r(count1,2)={decide_l};
                        list(ind)=count1;
                    end
                    if count==num
                        criteria1=0;
                    end
                end
                
                s_box=[label_r{:,2}];
                s_box=s_box([true,diff(sort(s_box))>0]);
                for j=1:length(s_box)
                    remove_label_box=[label_r{logical([label_r{:,2}]==s_box(j))}];
                    Label_record{label+j,1}=cluster_volum(remove_label_box,:)';
                    Label_record{label+j,2}=cluster_volum2(remove_label_box)';
                    Lable_mat(cluster_volum2(remove_label_box))=label+j;
                    Label_record{label+j,3}=label+j;
                end
                
                label=label+length(s_box);
            else
                Label_record{label+1,1}=cluster_volum';
                Label_record{label+1,2}=cluster_volum2';
                Label_record{label+1,3}=label+1;
                Lable_mat(cluster_volum2)=label+1;
                label=label+1;        
            end
        end    
        
        fb=[];
        criteria=[];
        decide_label_1=[];
        check_volum_find=[];
        
end

end

        




