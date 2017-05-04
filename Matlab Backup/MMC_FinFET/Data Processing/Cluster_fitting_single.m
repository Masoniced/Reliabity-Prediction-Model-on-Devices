function [p,data,final,curve]=Cluster_fitting_single(Data,ratio,fixa,fixb,fixtor)
syms a1 b1 tor1 a2 b2 tor2 t
CDF=1-((1+1/a1*(t/tor1)^b1)^(-a1));
PDF=diff(CDF,t);
LF=-log(PDF);
% s=xlsread('C:\D\Program Files\MATLAB\projeects\Clustering data processing\FinFET\125 TDDB FR09.xlsx');
% Data=s(:,2);
data=Data;
data=data(~isnan(data));
data=sort(data);
x0=[0.001,1,data(round(length(data)*0.63*ratio))];
%x0=[1,1,5];
alph=5;
fa1=diff(LF,a1);
fb1=diff(LF,b1);
ftor1=diff(LF,tor1);
% fa2=diff(LF,a2);
% fb2=diff(LF,b2);
% ftor2=diff(LF,tor2);
Lf=matlabFunction(LF);
G1=matlabFunction(fa1);
% G2=matlabFunction(fa2);
G3=matlabFunction(fb1);
% G4=matlabFunction(fb2);
G5=matlabFunction(ftor1);
% G6=matlabFunction(ftor2);
t_G1=0;
t_G2=0;
t_G3=0;
t_G4=0;
t_G5=0;
t_G6=0;
F=0;

if isempty(fixa)&&isempty(fixb)&&isempty(fixtor)
    for i=1:length(data)
        t_G1=G1(a1,b1,data(i),tor1)+t_G1;
%         t_G2=G2(a1,a2,b1,b2,data(i),tor1,tor2)+t_G2;
        t_G3=G3(a1,b1,data(i),tor1)+t_G3;
%         t_G4=G4(a1,a2,b1,b2,data(i),tor1,tor2)+t_G4;
        t_G5=G5(a1,b1,data(i),tor1)+t_G5;
%         t_G6=G6(a1,a2,b1,b2,data(i),tor1,tor2)+t_G6;        
        F=Lf(a1,b1,data(i),tor1)+F;
    end
    t_G1=matlabFunction(t_G1);
%     t_G2=matlabFunction(t_G2);
    t_G3=matlabFunction(t_G3);
%     t_G4=matlabFunction(t_G4);
    t_G5=matlabFunction(t_G5);
%     t_G6=matlabFunction(t_G6);    
    F=matlabFunction(F);

    temp_x=x0;
    criteria=1e-10;
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
            grad=[t_G1(temp_x(1),temp_x(2),temp_x(3)),t_G3(temp_x(1),temp_x(2),temp_x(3)),t_G5(temp_x(1),temp_x(2),temp_x(3))];
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
    high_bond=max(data)*3;
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
        pre_test=CDF(final(1),final(2),i,final(3));
%         if pre_test>=low_limit&&pre_test<=high_limit
        test=[test,pre_test];
%         end
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
elseif ~isempty(fixb)&&isempty(fixa)&&isempty(fixtor)
    x0=[0.001,1,data(round(length(data)*0.63*ratio))];
    for i=1:length(data)
        t_G1=G1(a1,a2,fixb(1),fixb(2),data(i),tor1,tor2)+t_G1;
        t_G2=G2(a1,a2,fixb(1),fixb(2),data(i),tor1,tor2)+t_G2;
        t_G5=G5(a1,a2,fixb(1),fixb(2),data(i),tor1,tor2)+t_G5;
        t_G6=G6(a1,a2,fixb(1),fixb(2),data(i),tor1,tor2)+t_G6;        
        F=Lf(a1,a2,fixb(1),fixb(2),data(i),tor1,tor2)+F;
    end
    t_G1=matlabFunction(t_G1);
    t_G2=matlabFunction(t_G2);
    t_G5=matlabFunction(t_G5);
    t_G6=matlabFunction(t_G6);    
    F=matlabFunction(F);

    temp_x=x0;
    criteria=1e-8;
    error=1;
    c=1;
    F1=F(temp_x(1),temp_x(2),temp_x(3),temp_x(4));
    count=0;
    list=[];
    while error>criteria
        count=count+1;
        list(count,1:4)=temp_x;
        list(count,5)=F1;
        if c==1
            grad=[t_G1(temp_x(1),temp_x(2),temp_x(3),temp_x(4)),t_G2(temp_x(1),temp_x(2),temp_x(3),temp_x(4)),t_G5(temp_x(1),temp_x(2),temp_x(3),temp_x(4)),t_G6(temp_x(1),temp_x(2),temp_x(3),temp_x(4))];
            p=grad./sqrt(grad(1)^2+grad(2)^2+grad(3)^2+grad(4)^2);
        end
        next_x=temp_x-alph*p;
        if next_x(1)<0
            F2=F1+1;
        elseif next_x(2)<0
            F2=F1+1;
        elseif next_x(3)<0
            F2=F1+1;
        elseif next_x(4)<0
            F2=F1+1;
        else
            F2=F(next_x(1),next_x(2),next_x(3),next_x(4));
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

    final=[temp_x(1),temp_x(2),temp_x(3),temp_x(4)];
    CDF=matlabFunction(CDF);
    % for i=1:length(data)
    %     prob(i)=CDF(final(1),final(2),data(i),final(3));
    %     y(i)=log(-log(1-prob(i)));
    % end
    test=[];
    ccc=0;
    count1=0;
    low_bond=min(data)/4;
    high_bond=max(data)*10;
    low_limit=1-exp(-exp(-8));
    high_limit=1-exp(-exp(3.5));
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
        pre_test=CDF(final(1),final(2),fixb(1),fixb(2),i,final(3),final(4));
%         if pre_test>=low_limit&&pre_test<=high_limit
        test=[test,pre_test];
%         end
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    curve=[j(ind);test(ind)];
    
    n=1:length(data);
    p=(n-0.3)/(length(data)+0.4);
    p=log(-log(1-p));
    
elseif ~isempty(fixa)&&isempty(fixb)&&isempty(fixtor)
    x0=[1,1,data(round(length(data)*0.63*ratio)),data(round(length(data)*0.63*(1-ratio)+length(data)*ratio))];
    for i=1:length(data)
        t_G3=G3(fixa(1),fixa(2),b1,b2,data(i),tor1,tor2)+t_G3;
        t_G4=G4(fixa(1),fixa(2),b1,b2,data(i),tor1,tor2)+t_G4;
        t_G5=G5(fixa(1),fixa(2),b1,b2,data(i),tor1,tor2)+t_G5;
        t_G6=G6(fixa(1),fixa(2),b1,b2,data(i),tor1,tor2)+t_G6;        
        F=Lf(fixa(1),fixa(2),b1,b2,data(i),tor1,tor2)+F;
    end
    t_G3=matlabFunction(t_G3);
    t_G4=matlabFunction(t_G4);
    t_G5=matlabFunction(t_G5);
    t_G6=matlabFunction(t_G6);    
    F=matlabFunction(F);

    temp_x=x0;
    criteria=1e-8;
    error=1;
    c=1;
    F1=F(temp_x(1),temp_x(2),temp_x(3),temp_x(4));
    count=0;
    list=[];
    while error>criteria
        count=count+1;
        list(count,1:4)=temp_x;
        list(count,5)=F1;
        if c==1
            grad=[t_G3(temp_x(1),temp_x(2),temp_x(3),temp_x(4)),t_G4(temp_x(1),temp_x(2),temp_x(3),temp_x(4)),t_G5(temp_x(1),temp_x(2),temp_x(3),temp_x(4)),t_G6(temp_x(1),temp_x(2),temp_x(3),temp_x(4))];
            p=grad./sqrt(grad(1)^2+grad(2)^2+grad(3)^2+grad(4)^2);
        end
        next_x=temp_x-alph*p;
        if next_x(1)<0
            F2=F1+1;
        elseif next_x(2)<0
            F2=F1+1;
        elseif next_x(3)<0
            F2=F1+1;
        elseif next_x(4)<0
            F2=F1+1;
        else
            F2=F(next_x(1),next_x(2),next_x(3),next_x(4));
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

    final=[temp_x(1),temp_x(2),temp_x(3),temp_x(4)];
    CDF=matlabFunction(CDF);
    % for i=1:length(data)
    %     prob(i)=CDF(final(1),final(2),data(i),final(3));
    %     y(i)=log(-log(1-prob(i)));
    % end
    test=[];
    ccc=0;
    count1=0;
    low_bond=min(data)/4;
    high_bond=max(data)*10;
    low_limit=1-exp(-exp(-4));
    high_limit=1-exp(-exp(3.5));
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
        pre_test=CDF(fixa(1),fixa(2),final(1),final(2),i,final(3),final(4));
%         if pre_test>=low_limit&&pre_test<=high_limit
        test=[test,pre_test];
%         end
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    curve=[j(ind);test(ind)];
    
    n=1:length(data);
    p=(n-0.3)/(length(data)+0.4);
    p=log(-log(1-p));
elseif ~isempty(fixa)&&(~isempty(fixb))&&isempty(fixtor)
    x0=data(round(length(data)*0.63*ratio));
    for i=1:length(data)
%         t_G1=G1(a1,b1,data(i),tor1)+t_G1;
%         t_G2=G2(a1,a2,b1,b2,data(i),tor1,tor2)+t_G2;
%         t_G3=G3(a1,b1,data(i),tor1)+t_G3;
%         t_G4=G4(a1,a2,b1,b2,data(i),tor1,tor2)+t_G4;
        t_G5=G5(fixa,fixb,data(i),tor1)+t_G5;
%         t_G6=G6(a1,a2,b1,b2,data(i),tor1,tor2)+t_G6;        
        F=Lf(fixa,fixb,data(i),tor1)+F;
    end
%     t_G1=matlabFunction(t_G1);
%     t_G2=matlabFunction(t_G2);
%     t_G3=matlabFunction(t_G3);
%     t_G4=matlabFunction(t_G4);
    t_G5=matlabFunction(t_G5);
%     t_G6=matlabFunction(t_G6);    
    F=matlabFunction(F);

    temp_x=x0;
    criteria=1e-10;
    error=1;
    c=1;
    F1=F(temp_x);
    count=0;
    list=[];
    while error>criteria
        count=count+1;
        list(count,1)=temp_x;
        list(count,2)=F1;
        if c==1
            grad=t_G5(temp_x(1));
            p=grad./sqrt(grad^2);
        end
        next_x=temp_x-alph*p;
        if next_x(1)<0
            F2=F1+1;
        else
            F2=F(next_x(1));
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

    final=temp_x;
    CDF=matlabFunction(CDF);
    % for i=1:length(data)
    %     prob(i)=CDF(final(1),final(2),data(i),final(3));
    %     y(i)=log(-log(1-prob(i)));
    % end
    test=[];
    ccc=0;
    count1=0;
    low_bond=min(data)/2;
    high_bond=max(data)*3;
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
        pre_test=CDF(fixa,fixb,i,final);
%         if pre_test>=low_limit&&pre_test<=high_limit
        test=[test,pre_test];
%         end
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    curve=[j(ind);test(ind)];
    
    n=1:length(data);
    p=(n-0.3)/(length(data)+0.4);
    p=log(-log(1-p));
else
    CDF=matlabFunction(CDF);
    test=[];
    ccc=0;
    count1=0;
    low_bond=min(data)/4;
    high_bond=max(data)*10;
    low_limit=1-exp(-exp(-4));
    high_limit=1-exp(-exp(3.5));
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
        pre_test=CDF(fixa(1),fixa(2),fixb(1),fixb(2),i,fixtor(1),fixtor(2));
%         if pre_test>=low_limit&&pre_test<=high_limit
        test=[test,pre_test];
%         end
    end
    test=log(-log(1-test));
    ind=~isinf(test);
    curve=[j(ind);test(ind)];
    
    n=1:length(data);
    p=(n-0.3)/(length(data)+0.4);
    p=log(-log(1-p));     
    final=[0,0];
end

% figure
% hold on
% plot(curve(1,:),curve(2,:))
% plot(data,p,'*')
% hold off
end