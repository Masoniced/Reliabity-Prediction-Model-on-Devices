syms a1 b1 x a2 b2 y t
CDF=1-((1+1/a1*(t/x)^b1)^(-a1))*((1+1/a2*(t/y)^b2)^(-a2));


s1=xlsread('C:\D\Program Files\MATLAB\projeects\Clustering data processing\FinFET\125 TDDB FR2729.xlsx');
s2=xlsread('C:\D\Program Files\MATLAB\projeects\Clustering data processing\FinFET\125 TDDB FR0709.xlsx');
s3=xlsread('C:\D\Program Files\MATLAB\projeects\Clustering data processing\FinFET\125 TDDB FR37.xlsx');

num=6;
Data={};
data={};
Data{1}=s1(:,1);
Data{2}=s1(:,9);
Data{3}=s2(:,5);
Data{4}=s2(:,9);
Data{5}=s3(:,1)*log(5);
Data{6}=s3(:,5)*log(5);

G=[];
LossF=[];
p_list={};
x0_final=[0.001,0.001,1,1];
ratio=1/3;    

G=sym(zeros(1,2*num+4));
for i=1:num
    [t_G1,t_G2,t_G3,t_G4,t_G5,t_G6,F,x0,data_s,p]=Cluster_combination([Data{i}],ratio,i);
    if i==1
        G(1)=t_G1;
        G(2)=t_G2;
        G(3)=t_G3;
        G(4)=t_G4;
        G(5)=t_G5;
        G(6)=t_G6;
        LossF=F;
    else
        G(1)=G(1)+t_G1;
        G(2)=G(2)+t_G2;
        G(3)=G(3)+t_G3;
        G(4)=G(4)+t_G4;
        G(2*i+3)=t_G5;
        G(2*i+4)=t_G6;
        LossF=LossF+F;
    end
    x0_final=[x0_final,x0];
    p_list(i)={p};
    data(i)={data_s};
end

G=matlabFunction(G);
LossF=matlabFunction(LossF);

x0_final=num2cell(x0_final);

temp_x=x0_final;
criteria=1e-10;
error=1;
c=1;
alph=5;
F1=LossF(temp_x{:});
count=0;
list=[];
while error>criteria
    count=count+1;
    list(count,1:2*num+4)=[temp_x{:}];
    list(count,2*num+5)=F1;
    if c==1
        grad=G(temp_x{:});        
        p=grad./sqrt(sum(grad.^2));
    end
    next_x=[temp_x{:}]-alph*p;
    next_x=num2cell(next_x);
    if next_x{1}<0
        F2=F1+1;
    elseif next_x{2}<0
        F2=F1+1;
    elseif next_x{3}<0
        F2=F1+1;
    else
        F2=LossF(next_x{:});
    end
    if F2<F1
        temp_x=next_x;
        error=abs((F2-F1)/F1);
        F1=F2;
        c=1;
%             alph=5;
    else
        alph=alph/2;
        c=0;
    end
    if alph<criteria
        error=0;
    end
end

final=[temp_x{:}];
vector=cell(1,num);
for i=1:num
    X_factor=subs(CDF,a1,final(1));
    X_factor=subs(X_factor,a2,final(2));
    X_factor=subs(X_factor,b1,final(3));
    X_factor=subs(X_factor,b2,final(4));
    X_factor=subs(X_factor,x,final(2*i+3));
    X_factor=subs(X_factor,y,final(2*i+4));
    X_factor=matlabFunction(X_factor);

    test=[];
    ccc=0;
    count1=0;
    low_bond=min([Data{i}])/20;
    high_bond=max([Data{i}])*200;
    low_limit=1-exp(-exp(-4));
    high_limit=1-exp(-exp(3.5));
    leng=(max([Data{i}])-min([Data{i}]));
    j=[];
    while ccc==0
        j=[j,low_bond+10*(count1)*leng/200];
        count1=count1+1;
        if j(end)>high_bond
            j(end)=[];
            ccc=1;
        end
    end
    for u=j
        pre_test=X_factor(u);
%         if pre_test>=low_limit&&pre_test<=high_limit
        test=[test,pre_test];
%         end
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    vector(i)={[j(ind);test(ind)]};
end

figure
hold on
plot(vector{1}(1,:),vector{1}(2,:),'r')
plot([data{1}],[p_list{1}],'*')
plot(vector{2}(1,:),vector{2}(2,:),'g')
plot([data{2}],[p_list{2}],'+')
plot(vector{3}(1,:),vector{3}(2,:),'b')
plot([data{3}],[p_list{3}],'*')
plot(vector{4}(1,:),vector{4}(2,:),'c')
plot([data{4}],[p_list{4}],'+')
plot(vector{5}(1,:),vector{5}(2,:),'y')
plot([data{5}],[p_list{5}],'*')
plot(vector{6}(1,:),vector{6}(2,:),'m')
plot([data{6}],[p_list{6}],'+')
hold off
    
    