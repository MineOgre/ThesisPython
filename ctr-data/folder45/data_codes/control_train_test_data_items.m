function [ a ] = control_train_test_data_items( v_train, v_test, items )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

R_tr = gen_R_from_items(v_train);
R_ts = gen_R_from_items(v_test);
R = gen_R_from_items(items);
R_sum = R_tr + R_ts;
a = isequal(R_sum,R);

end

