function [u,data,final,curve]=Least_fitting_square_clustering(Data,ratio,fixa,fixb)
syms a1 b1 tor1 a2 b2 tor2 x
CDF=1-((1+1/a1*(x/tor1)^b1)^(-a1))*((1+1/a2*(x/tor2)^b2)^(-a2));

data=Data;
data=data(~isnan(data));
data=sort(data);

alph=5;

F=matlabFunction(CDF);
error_Function=0;

n=1:length(data);
p=(n-0.3)/(length(data)+0.4);
u=log(-log(1-p));


if isempty(fixa)&&isempty(fixb)
    input_vec=[a1,a2,b1,b2,tor1,tor2];
    n_input=[a1,a2,b1,b2,tor1,tor2];
    x0=[0.001,0.001,1,1,data(round(length(data)*0.63*ratio)),data(round(length(data)*0.63*(1-ratio)+length(data)*ratio))];
    CDF_X=CDF;
elseif isempty(fixa)&&(~isempty(fixb))
    input_vec=[a1,a2,fixb(1),fixb(2),tor1,tor2];
    n_input=[a1,a2,tor1,tor2];
    x0=[0.001,0.001,data(round(length(data)*0.63*ratio)),data(round(length(data)*0.63*(1-ratio)+length(data)*ratio))];
    CDF_X=subs(CDF,b1,fixb(1));
    CDF_X=subs(CDF_X,b2,fixb(2));
elseif (~isempty(fixa))&&isempty(fixb)
    input_vec=[fixa(1),fixa(2),b1,b2,tor1,tor2];
    n_input=[b1,b2,tor1,tor2];
    x0=[1,1,data(round(length(data)*0.63*ratio)),data(round(length(data)*0.63*(1-ratio)+length(data)*ratio))];
    CDF_X=subs(CDF,a1,fixa(1));
    CDF_X=subs(CDF_X,a2,fixa(2));
elseif (~isempty(fixa))&&(~isempty(fixb))
    input_vec=[fixa(1),fixa(2),fixb(1),fixb(2),tor1,tor2];
    n_input=[tor1,tor2];
    x0=[data(round(length(data)*0.63*ratio)),data(round(length(data)*0.63*(1-ratio)+length(data)*ratio))];
    CDF_X=subs(CDF,b1,fixb(1));
    CDF_X=subs(CDF_X,b2,fixb(2));    
    CDF_X=subs(CDF_X,a1,fixa(1));
    CDF_X=subs(CDF_X,a2,fixa(2));    
end
    

for i=1:length(data)
    vec=num2cell([input_vec,data(i)]);
    error_Function=error_Function+(F(vec{:})-p(i))^2;
end

G=[];
for i=1:length(n_input)
    G=[G,diff(error_Function,n_input(i))];
end

G=matlabFunction(G);
error_Function=matlabFunction(error_Function);

temp_x=num2cell(x0);
criteria=1e-8;
error=1;
c=1;
F1=error_Function(temp_x{:});
count=0;
list=[];

while error>criteria
    list=[list,F1];
    count=count+1;
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
    else
        F2=error_Function(next_x{:});
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
X_factor=CDF_X;
for i=1:length(final);
    X_factor=subs(X_factor,n_input(i),final(i));
end

X_factor=matlabFunction(X_factor);

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
    pre_test=X_factor(i);
%         if pre_test>=low_limit&&pre_test<=high_limit
    test=[test,pre_test];
%         end
end
test=log(-log(1-test));
ind=~isinf(test);
curve=[j(ind);test(ind)];

end

