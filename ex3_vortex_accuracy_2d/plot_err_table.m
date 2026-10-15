%Plot table 
clear all; close all;

format short e

err_load =load('err_table.dat'); 


n_array = err_load(:,1)';
L1Error = err_load(:,2);
L1Order = err_load(:,3);
LiError = err_load(:,4);
LiOrder = err_load(:,5);



%Format table
fractions = arrayfun(@(x) sprintf('1/%d', x), n_array, 'UniformOutput', false);
T = table(fractions', L1Error(:),L1Order(:), LiError(:),LiOrder(:),...
    'VariableNames',{'dx','L1 Error','L1 Order', 'Li Error', 'Li Order'});
disp(T)