function [new_Probability]=Probability_vector_process(process_Probability,ind)
% pick_y==1 Generation
% pick_y==2 Recombination
% pick_y==3 Hopping
new_Probability=process_Probability;
new_Probability(ind,:)=[];
vector=process_Probability(:,1);
L=length(vector);
if L==1
    new_Probability=[];
else
    if ind==1
        vector=vector-vector(ind);
        vector(ind)=[];
        vector=vector/vector(end);
        new_Probability(:,1)=vector;
    elseif ind==L
        vector(ind)=[];
        vector=vector/vector(ind-1);
        new_Probability(:,1)=vector;
    else
        vector(ind+1:end)=vector(ind+1:end)-(vector(ind)-vector(ind-1));
        vector(ind)=[];
        vector=vector/vector(end);
        new_Probability(:,1)=vector;
    end
end
        