ss=Shape_record.plane_view(:,:,number1);
ss=ss*(5.43)/4*sqrt(2);
ss1=zeros(size(ss));
for i=1:544
    ss1(i,:)=ss(545-i,:);
end
