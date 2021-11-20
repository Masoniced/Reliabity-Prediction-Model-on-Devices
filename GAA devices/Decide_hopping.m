function [destination]=Decide_hopping(p)
global M n1 n2 n3
%tic
%% Decide Hopping Position
criteria=1;
test_p=p;
while criteria==1
    [k]=sample(test_p,n1,n2,n3-2,0);
    ind=k(M(k)==0);
    if isempty(ind)
        s=ceil(rand(1,1)*length(k));
        test_p=k(s);
    else
        f=ceil(rand(1,1)*length(ind));
        destination=ind(f);
        criteria=0;
    end
end
%toc
end
    