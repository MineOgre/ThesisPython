function [ R ] = gen_R_from_items( items )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
num_v = size(items,1);
num_u= max(max(items(:,2:end)))+1;

R = zeros(num_u, num_v);
for i=1:num_v 
    for j = 2:(items(i,1)+1)
        R(items(i,j)+1,i)=1;
    end
end

end

