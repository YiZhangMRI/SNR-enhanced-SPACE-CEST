function [echo_amp,alpha_array] = func_simu_pss(T1,T2,tesp,ETL,terminal_FA)
% simulate pseudo-steady state signal evolution

%% input
% T1 = 1000/1000;
% T2 = 100/1000;
% tesp = 3.46/1000;
% ETL = 140;
% terminal_FA = 120;


%% flip angle look-up table
num = 8;
FAs = [ 148.9, 121.8, 118.8, 120.4, 120.1, 119.7, 120.1, 120.0 ];     % 120 degrees
  


%%
assert(terminal_FA==120,'only support 120 degrees')
alpha_array = terminal_FA*ones(1, ETL);
for k = 1:min([num, ETL])
    alpha_array(k) = FAs(k);
end
phi_array = zeros(size(alpha_array));
echo_amp = func_forward_EPG(T1, T2, tesp, alpha_array, phi_array);
end