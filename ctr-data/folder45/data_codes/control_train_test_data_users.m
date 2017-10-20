function [ a ] = control_train_test_data_users( u_train, u_test, users, v_num )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

R_tr = gen_R_from_users(u_train, v_num);
R_ts = gen_R_from_users(u_test, v_num);
R = gen_R_from_users(users, v_num);
R_sum = R_tr + R_ts;
a = isequal(R_sum,R);

end

