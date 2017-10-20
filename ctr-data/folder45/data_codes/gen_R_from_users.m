function [ R ] = gen_R_from_users( users, v_num )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u_num = size(users,1);

R = zeros(u_num, v_num);
for i=1:u_num 
    for j = 2:(users(i,1)+1)
        %if (users(i,j)~=0)
            R(i,users(i,j)+1)=1;
        %end
    end
end

end

