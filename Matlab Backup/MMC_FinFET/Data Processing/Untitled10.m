%% Fourier analysis
Data=load('C:\D\Program Files\MATLAB\projeects\Clustering data processing\Fourier\current.mat');
D=Data.res;
list=cell(size(D));
example=cell(1,length(D));
for i=1:numel(D)
    test=D{i};
    [ix,iy]=ind2sub(size(D),i);
    if isempty(test)
        list(ix,iy)={[]};
    else
        
        factor=3.5;
        L=length(test(:,1));
        I=test(:,2);
        S_L=2^nextpow2(L);
        F_t=fft(I,S_L);
        dt=1/(abs(test(1,1)-test(2,1)));
        f=dt/2*linspace(0,1,S_L/2+1);
        u=f;
        u(1)=[];
        ind=round((S_L/2+1)/factor);
        f(1:ind)=[];
        f(end-ind:end)=[];
        f=log(f);
        p=F_t.*conj(F_t)/S_L;
        f_t=log(p(ind+1:S_L/2-ind))';
        A=zeros(length(f_t),2);
        A(:,1)=f';
        A(:,2)=1;
        b=f_t';
        [x,~]=lsqr(A,b,1e-12,500);
        list(ix,iy)={x(1)};            
    end
    if ix==1
        example(iy)={[u;p(2:S_L/2+1)']};
    end
end

fs=6;
low_list1=[list{1,:}];
low_list1=sort(low_list1);
ind1=round(length(low_list1)/fs);
low_list1(1:ind1)=[];
low_list1(end-ind1:end)=[];

low_list2=[list{4,:}];
low_list2=sort(low_list2);
ind2=round(length(low_list2)/fs);
low_list2(1:ind2)=[];
low_list2(end-ind2:end)=[];

high_list1=[list{2,:}];
high_list1=sort(high_list1);
ind3=round(length(high_list1)/fs);
high_list1(1:ind3)=[];
high_list1(end-ind3:end)=[];

high_list2=[list{3,:}];
high_list2=sort(high_list2);
ind4=round(length(high_list2)/fs);
high_list2(1:ind4)=[];
high_list2(end-ind4:end)=[];

high_list3=[list{5,:}];
high_list3=sort(high_list3);
ind5=round(length(high_list3)/fs);
high_list3(1:ind5)=[];
high_list3(end-ind5:end)=[];

mean_list=[mean(-low_list1),mean(-low_list2),mean(-high_list1-0.55),mean(-high_list2),mean(-high_list3)];
std_list=[var(-low_list1),var(-low_list2),var(-high_list1),var(-high_list2),var(-high_list3)];

% figure 
% hold on
% errorbar(2,mean_list(1),std_list(1))
% errorbar(2.05,mean_list(2),std_list(2))
% errorbar(2.1,mean_list(3),std_list(3))
% errorbar(2.15,mean_list(4),std_list(4))
% errorbar(2.2,mean_list(5),std_list(5))
% % scatter(abs(low_list1),ones(length(low_list1),1)*2)
% % scatter(abs(low_list2),ones(length(low_list2),1)*2.05)
% % scatter(abs(high_list1+0.55),ones(length(high_list1),1)*2.1)
% % scatter(abs(high_list2),ones(length(high_list2),1)*2.15)
% % scatter(abs(high_list3),ones(length(high_list3),1)*2.2)
% hold off

figure
hold on
for i=1:length(example)
    plot(example{i}(1,:),example{i}(2,:))
end
hold off