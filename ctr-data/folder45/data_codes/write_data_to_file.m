function [ ] = write_data_to_file( data, file_l )
%WRITE_DATA_TO_FILE Summary of this function goes here
%   Detailed explanation goes here

[u,v] = size(data);

for i=1:u
    num = data(i,1);
    tmp = data(i,1:(num+1));
    if (i==1)
        dlmwrite(file_l,tmp, 'delimiter',' ');
    else
        dlmwrite(file_l,tmp, 'delimiter',' ','-append');
    end
    clear tmp;
end

end

