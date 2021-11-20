%% Random walk path generation
function [R,Check_value] = RPG_r(criteria,PIT_p,n1,n2,n3,t)

R = zeros(n1,n2,n3);
if criteria == 1
    fprintf('Criteria error');
    return
end
Check_value = [];
rng('shuffle')
x = round(rand(1,1)*n1);
if x<1
    x = 1;
end
initial_p=[x,1,1];
R(x,1,1)=1;
Check_value=[x,1,1];
label_value=cell(1,n3-t);
label_value(1)={[x,1,1]};
direction=[-1,1,0;1,1,0;0,1,0];
check=0;
check1=0;
for i=1:n3-t
    start_p=initial_p;
    rng('shuffle')
    if i==1
        while check==0
            rnumber=round(3*rand(1,1));
            if rnumber<1
                rnumber=1;
            end
            if start_p(2)==1
                residul=[start_p(1),start_p(2),start_p(3)];
            end
            next_p=start_p+direction(rnumber,:);
            if next_p(1)==0
                next_p(1)=1;
            end
            if next_p(1)==n1+1
                next_p(1)=n1;
            end
            R(next_p(1),next_p(2),next_p(3))=1;
            Check_value=[Check_value;next_p(1),next_p(2),next_p(3)];
            label_value(i)={[label_value{i};next_p(1),next_p(2),next_p(3)]};
            start_p=next_p;
            if next_p(2)==n2
                check=1;
                rnumber_2=round(3*rand(1,1));
                if rnumber_2<1
                    rnumber_2=1;
                end
                initial_p=residul+direction(rnumber_2,:);
                if initial_p(1)==0
                    initial_p(1)=1;
                end
                if initial_p(1)==n1+1
                    initial_p(1)=n1;
                end
                initial_p(2)=initial_p(2)-1;
                initial_p(3)=initial_p(3)+1;
                Check_value=[Check_value;initial_p(1),initial_p(2),initial_p(3)];
                label_value(i+1)={[label_value{i+1};initial_p(1),initial_p(2),initial_p(3)]};
            end
        end
    else
        reference=label_value{i-1};
        while check1==0;
            rnumber=round(3*rand(1,1));
            if rnumber<1
                rnumber=1;
            end
            if start_p(2)==1
                residul=[start_p(1),start_p(2),start_p(3)];
            end            
            o=reference(reference(:,2)==start_p(2)+1);
            if (abs(o-start_p(1)))==0
                next_p=start_p+direction(rnumber,:);
                if next_p(1)==0
                    next_p(1)=1;
                end
                if next_p(1)==n1+1
                    next_p(1)=n1;
                end
                R(next_p(1),next_p(2),next_p(3))=1;
                Check_value=[Check_value;next_p(1),next_p(2),next_p(3)];
                label_value(i)={[label_value{i};next_p(1),next_p(2),next_p(3)]};
                start_p=next_p;
            elseif abs((o-start_p(1)))==1
                if round(rand(1,1))==0
                    next_p=[start_p(1),start_p(2)+1,i];
                else
                    next_p=[o(1),start_p(2)+1,i];
                end
                if next_p(1)==0
                    next_p(1)=1;
                end
                if next_p(1)==n1+1
                    next_p(1)=n1;
                end
                R(next_p(1),next_p(2),next_p(3))=1;
                Check_value=[Check_value;next_p(1),next_p(2),next_p(3)];
                label_value(i)={[label_value{i};next_p(1),next_p(2),next_p(3)]};
                start_p=next_p;
            elseif abs((o-start_p(1)))==2
                if o>start_p(1)
                    next_p=[start_p(1)+1,start_p(2)+1,i];
                else
                    next_p=[start_p(1)-1,start_p(2)+1,i];     
                end
                if next_p(1)==0
                    next_p(1)=1;
                end
                if next_p(1)==n1+1
                    next_p(1)=n1;
                end
                R(next_p(1),next_p(2),next_p(3))=1;
                Check_value=[Check_value;next_p(1),next_p(2),next_p(3)];
                label_value(i)={[label_value{i};next_p(1),next_p(2),next_p(3)]};
                start_p=next_p;
            else
                fprintf('Path error')
                break
            end
            if next_p(2)==n2
                check1=1;
                rnumber_1=round(3*rand(1,1));
                if rnumber_1<1
                    rnumber_1=1;
                end
                initial_p=residul+direction(rnumber_1,:);
                if initial_p(1)==0
                    initial_p(1)=1;
                end
                if initial_p(1)==n1+1
                    initial_p(1)=n1;
                end
                initial_p(2)=initial_p(2)-1;
                initial_p(3)=initial_p(3)+1;
                if initial_p(3)<n3-t+1
                    Check_value=[Check_value;initial_p(1),initial_p(2),initial_p(3)];
                    label_value(i+1)={[label_value{i+1};initial_p(1),initial_p(2),initial_p(3)]};
                end
            end           
        end
        check1=0;
    end
end
 
rng('shuffle')
collect=[];
num_PIT=round(PIT_p*n1*n2*n3)-n1*n3;
if num_PIT<0
    fprintf('GB Error')
end
count=0;
ct=0;
while ct==0
    k=rand(1,1);
    ind = round(k*n1*n2*n3);
    if ind<1
        ind = ind+1;
    end
    if R(ind)==1
        continue
    else
       R(ind)=1;
       X(ind)=1;
       count=count+1;
    end
    if count==num_PIT
        ct=1;
    end
end
end