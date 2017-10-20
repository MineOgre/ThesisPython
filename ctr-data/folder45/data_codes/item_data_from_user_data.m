function [ items ] = item_data_from_user_data( users, u_num, v_num, size_v)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

items = zeros(v_num, size_v);
for i=1:u_num
    tmp = users(i,2:end)+1;
    i_num = users(i,1);
    for j=1:i_num  
        indx = items(tmp(j),1);
        items(tmp(j),indx+2) = i-1;
        items(tmp(j),1) = items(tmp(j),1)+1;
    end
end

end

