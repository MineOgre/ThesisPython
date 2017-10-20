clear;

%users = dlmread('./ctr-data/folder45/users_small.dat');
users = dlmread('../users.dat');
u_num = size(users,1);
v_num= max(max(users))+1;

P=10;
try 
    P;
catch
    disp('P value is missing!!!');
end

R = gen_R_from_users( users, v_num );
max_user_per_item = max(sum(R));
max_item_per_user = max(sum(R,2));
size_v = max_user_per_item+1;
size_u = max_item_per_user+1;

%%%%%
%%%%% Buras? users bilgisinden item matrisini ??karma k?sm?
items = item_data_from_user_data(users, u_num, v_num, size_v);
        
%dlmwrite('items_o.dat',items);
        
%%%%
%%%% sample training ve test set hazirlama from users
u_test_dim = max(users(:,1))-P;
u_train = zeros(u_num,P);
u_test = zeros(u_num,u_test_dim);

for i=1:u_num
    tmp = users(i,2:end);
    indx = users(i,1); 
    val = randsample(indx,P);
    u_train(i,1:P)=tmp(val);
    t = tmp;
    t(val) = [];
    u_test(i,:) = t;
end

u_train = [ones(u_num,1)*P u_train];   
u_test = [users(:,1)-P u_test];

%%%%
%%%% items sample training ve test set hazirlama 

R_train = gen_R_from_users( u_train, v_num );
train_max_user_per_item = max(sum(R_train));
train_max_item_per_user = max(sum(R_train,2));
train_size_v = train_max_user_per_item+1;
train_size_u = train_max_item_per_user+1;

test_size_v = max(sum(R-R_train))+1;

v_train = item_data_from_user_data(u_train, u_num, v_num, train_size_v);
v_test = item_data_from_user_data(u_test, u_num, v_num, test_size_v);

if (~control_train_test_data_users(u_train, u_test, users, v_num))
    disp('Hazirlanan training ve test datasi ana veriyle uyumlu degil-users')
end

if (~control_train_test_data_items(v_train, v_test, items))
    disp('Hazirlanan training ve test datasi ana veriyle uyumlu degil-items')
end
i=48;
%dlmwrite(sprintf('../../folder50/cf-train-%d-items.dat',P),v_train, 'delimiter','\t');
write_data_to_file(v_train,sprintf('../../folder%d/cf-train-%d-items.dat',i,P));
write_data_to_file(v_test,sprintf('../../folder%d/cf-test-%d-items.dat',i,P));
write_data_to_file(u_test,sprintf('../../folder%d/cf-test-%d-users.dat',i,P));
write_data_to_file(u_train,sprintf('../../folder%d/cf-train-%d-users.dat',i,P));
% dlmwrite(sprintf('../../folder50/cf-test-%d-items.dat',P),v_test, 'delimiter','\t');
% dlmwrite(sprintf('../../folder50/cf-test-%d-users.dat',P),u_test, 'delimiter','\t');
% dlmwrite(sprintf('../../folder50/cf-train-%d-users.dat',P),u_train, 'delimiter','\t');

clear i j i_num indx u_test_dim tmp t R_train val ;
        
        
        