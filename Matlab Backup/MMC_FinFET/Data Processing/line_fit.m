function [p]=line_fit(x,y)
syms a b
f=0;
for i=1:length(x)
    f=f+(a*x+b-y).^2;
end
fa=matlabFunction(diff(f,a));
fb=matlabFunction(diff(f,b));
F=matlabFunction(f);
error=1e-8;
r=1;
temp_p=[0,0];
F1=F(temp_p(1),temp_p(2));
alpha=5;
while r>error
    g=[fa(temp_p(1),temp_p(2)),fb(temp_p(1),temp_p(2))];
    next_p=temp_p-g*alpha;
    F2=F(next_p(1),next_p(2));
    if F2<F1
        temp_p=next_p;
        r=abs(F2-F1)/F1;
        F1=F2;
        alpha=5;
    else
        alpha=alpha/2;
    end
    if alpha<1e-6
        r=1e-9;
    end
end
p=temp_p;
end