function [p,data,final,curve]=Cluster_fitting(Data,fixa,fixb)
syms a b t tor
CDF=1-(1+1/a*(t/tor)^b)^(-a);
PDF=diff(CDF,t);
LF=-log(PDF);
% s=xlsread('C:\D\Program Files\MATLAB\projeects\Clustering data processing\FinFET\125 TDDB FR09.xlsx');
% Data=s(:,2);
data=Data;
data=data(~isnan(data));
data=sort(data);
x0=[0.001,1,data(round(length(data)*0.63))];
%x0=[1,1,5];
alph=5;
fa=diff(LF,a);
fb=diff(LF,b);
ftor=diff(LF,tor);
Lf=matlabFunction(LF);
G1=matlabFunction(fa);
G2=matlabFunction(fb);
G3=matlabFunction(ftor);
t_G1=0;
t_G2=0;
t_G3=0;
F=0;
if isempty(fixa)&&isempty(fixb)
    for i=1:length(data)
        t_G1=G1(a,b,data(i),tor)+t_G1;
        t_G2=G2(a,b,data(i),tor)+t_G2;
        t_G3=G3(a,b,data(i),tor)+t_G3;
        F=Lf(a,b,data(i),tor)+F;
    end
    t_G1=matlabFunction(t_G1);
    t_G2=matlabFunction(t_G2);
    t_G3=matlabFunction(t_G3);
    F=matlabFunction(F);

    temp_x=x0;
    criteria=1e-8;
    error=1;
    c=1;
    F1=F(temp_x(1),temp_x(2),temp_x(3));
    count=0;
    list=[];
    while error>criteria
        count=count+1;
        list(count,1:3)=temp_x;
        list(count,4)=F1;
        if c==1
            grad=[t_G1(temp_x(1),temp_x(2),temp_x(3)),t_G2(temp_x(1),temp_x(2),temp_x(3)),t_G3(temp_x(1),temp_x(2),temp_x(3))];
            p=grad./sqrt(grad(1)^2+grad(2)^2+grad(3)^2);
        end
        next_x=temp_x-alph*p;
        if next_x(1)<0
            F2=F1+1;
        elseif next_x(2)<0
            F2=F1+1;
        elseif next_x(3)<0
            F2=F1+1;
        else
            F2=F(next_x(1),next_x(2),next_x(3));
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

    final=[temp_x(1),temp_x(2),temp_x(3)];
    CDF=matlabFunction(CDF);
    % for i=1:length(data)
    %     prob(i)=CDF(final(1),final(2),data(i),final(3));
    %     y(i)=log(-log(1-prob(i)));
    % end
    test=[];
    ccc=0;
    count1=0;
    low_bond=min(data)/2;
    high_bond=max(data)*500;
    leng=(max(data)-min(data));
    j=[];
    while ccc==0
        j=[j,low_bond+10*(count1)*leng/500];
        count1=count1+1;
        if j(end)>high_bond
            j(end)=[];
            ccc=1;
        end
    end
    for i=j
        test=[test,CDF(final(1),final(2),i,final(3))];
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    curve=[j(ind);test(ind)];
    
    n=1:length(data);
    p=(n-0.3)/(length(data)+0.4);
    p=log(-log(1-p));

%     parmhat = wblfit(data);
%     p=wblcdf(data,parmhat(1),parmhat(2));
%     p=log(-log(1-p));
elseif ~isempty(fixb)
    x0=[0.001,data(round(length(data)*0.63))];
    for i=1:length(data)
        t_G1=G1(a,fixb,data(i),tor)+t_G1;
        t_G2=G2(a,fixb,data(i),tor)+t_G2;
        t_G3=G3(a,fixb,data(i),tor)+t_G3;
        F=Lf(a,fixb,data(i),tor)+F;
    end
    t_G1=matlabFunction(t_G1);
    t_G2=matlabFunction(t_G2);
    t_G3=matlabFunction(t_G3);
    F=matlabFunction(F);

    temp_x=x0;
    criteria=1e-8;
    error=1;
    c=1;
    F1=F(temp_x(1),temp_x(2));
    count=0;
    list=[];
    while error>criteria
        count=count+1;
        list(count,1:3)=temp_x;
        list(count,4)=F1;
        if c==1
            grad=[t_G1(temp_x(1),temp_x(2)),t_G3(temp_x(1),temp_x(2))];
            p=grad./sqrt(grad(1)^2+grad(2)^2);
        end
        next_x=temp_x-alph*p;
        if next_x(1)<0
            F2=F1+1;
        elseif next_x(2)<0
            F2=F1+1;
        else
            F2=F(next_x(1),next_x(2));
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

    final=[temp_x(1),temp_x(2)];
    CDF=matlabFunction(CDF);
    % for i=1:length(data)
    %     prob(i)=CDF(final(1),final(2),data(i),final(3));
    %     y(i)=log(-log(1-prob(i)));
    % end
    test=[];
    ccc=0;
    count1=0;
    low_bond=min(data)/2;
    high_bond=max(data)*500;
    leng=(max(data)-min(data));
    j=[];
    while ccc==0
        j=[j,low_bond+10*(count1)*leng/500];
        count1=count1+1;
        if j(end)>high_bond
            j(end)=[];
            ccc=1;
        end
    end
    for i=j
        test=[test,CDF(final(1),fixb,i,final(2))];
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    curve=[j(ind);test(ind)];
    
    n=1:length(data);
    p=(n-0.3)/(length(data)+0.4);
    p=log(-log(1-p));
elseif ~isempty(fixa)
        x0=[1,data(round(length(data)*0.63))];
    for i=1:length(data)
        t_G1=G1(fixa,b,data(i),tor)+t_G1;
        t_G2=G2(fixa,b,data(i),tor)+t_G2;
        t_G3=G3(fixa,b,data(i),tor)+t_G3;
        F=Lf(fixa,b,data(i),tor)+F;
    end
    t_G1=matlabFunction(t_G1);
    t_G2=matlabFunction(t_G2);
    t_G3=matlabFunction(t_G3);
    F=matlabFunction(F);

    temp_x=x0;
    criteria=1e-8;
    error=1;
    c=1;
    F1=F(temp_x(1),temp_x(2));
    count=0;
    list=[];
    while error>criteria
        count=count+1;        
        list(count,1:2)=temp_x;
        list(count,4)=F1;
        if c==1
            grad=[t_G2(temp_x(1),temp_x(2)),t_G3(temp_x(1),temp_x(2))];
            p=grad./sqrt(grad(1)^2+grad(2)^2);
        end
        next_x=temp_x-alph*p;
        if next_x(1)<0
            F2=F1+1;
        elseif next_x(2)<0
            F2=F1+1;
        else
            F2=F(next_x(1),next_x(2));
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

    final=[temp_x(1),temp_x(2)];
    CDF=matlabFunction(CDF);
    % for i=1:length(data)
    %     prob(i)=CDF(final(1),final(2),data(i),final(3));
    %     y(i)=log(-log(1-prob(i)));
    % end
    test=[];
    ccc=0;
    count1=0;
    low_bond=min(data)/2;
    high_bond=max(data)*500;
    leng=(max(data)-min(data));
    j=[];
    while ccc==0
        j=[j,low_bond+10*(count1)*leng/500];
        count1=count1+1;
        if j(end)>high_bond
            j(end)=[];
            ccc=1;
        end
    end
    for i=j
        test=[test,CDF(fixa,final(1),i,final(2))];
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    curve=[j(ind);test(ind)];
    
    n=1:length(data);
    p=(n-0.3)/(length(data)+0.4);
    p=log(-log(1-p));

%     parmhat = wblfit(data);
%     p=wblcdf(data,parmhat(1),parmhat(2));
%     p=log(-log(1-p));
end

% figure
% hold on
% plot(curve(1,:),curve(2,:))
% plot(data,p,'*')
% hold off
end



        
    

