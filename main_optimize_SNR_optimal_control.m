clear
close all

%% input
T1 = 1000/1000;
T2 = 100/1000;
tesp = 3.42/1e3;
nominal_turbo_factor = 140;
num_discard_echoes = 0;
num_total_echoes = nominal_turbo_factor + num_discard_echoes;
alpha_max = 180; % in degree
alpha_min = 0;

% SNR+Res
SNR_obj_weight = 1;
res_obj_weight = 2e-4;

% % Res-only
% SNR_obj_weight = 0;
% res_obj_weight = 1;



klim = 25;
num_channel = 1;
num_spatial_points = 2;
B1map = (1/num_channel) * ones(num_spatial_points, num_channel);


%%
tmpStruct = load('atten_rawdata_5subj_avg.mat');
atten = tmpStruct.atten;
atten = [zeros(1,num_discard_echoes), atten];
atten = atten ./ max(atten);
atten = atten(:).';

lb_oneCH = [90,alpha_min*ones(1, num_total_echoes)];
ub_oneCH = [90,alpha_max*ones(1, num_total_echoes)];
x0_oneCH = [90,10.*ones(1, num_total_echoes)];



%% use the notation of optimal control EPG
lb = repmat(lb_oneCH, 1, num_channel);
ub = repmat(ub_oneCH, 1, num_channel);
x0 = repmat(x0_oneCH, 1, num_channel);

lb = [lb, zeros(size(lb))]; % real, imag
ub = [ub, zeros(size(ub))];
x0 = [x0, zeros(size(x0))];

lb = deg2rad(lb);
ub = deg2rad(ub);
x0 = deg2rad(x0);
ESP = tesp;
W = [0,0,atten];
obj_fun = @(alpha_array)SNR_res_obj(alpha_array, T1,T2,ESP,W,SNR_obj_weight,res_obj_weight,B1map,klim);


%% do optimization
opt = optimoptions('fmincon', 'MaxIterations', 1000, 'Display', 'iter', ...
                   'SpecifyObjectiveGradient',true, 'UseParallel',true, ...
                   'GradConstr','off','DerivativeCheck','off');
opt = optimoptions(opt,'PlotFcns',@optimplotfval);               
[alpha_arr,fval,exitflag,output] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], opt);

    
%% show results
[~, ~,FF] = obj_fun(alpha_arr);
FF_1st_pixel = FF(:,:,1);
NN = size(FF_1st_pixel,1) / 2;
sigs = FF_1st_pixel(2,3:end) + 1i*FF_1st_pixel(NN+2,3:end);
sigs = abs(sigs);
% figure;
% plot(sigs);


alphas_all_channel = rad2deg(alpha_arr(1:numel(alpha_arr)/2));
alphas_1st_channel = alphas_all_channel(1:numel(alphas_all_channel)/num_channel);
refocus_alphas = alphas_1st_channel(2:end);


echo_amp = func_simu_pss(T1,T2,tesp,nominal_turbo_factor,120);
echo_times = tesp.*(1:num_total_echoes);
signal = func_forward_EPG(T1, T2, tesp, refocus_alphas);
figure;
subplot(121)
plot(exp(-echo_times./T2));
hold on;
plot(signal);
plot(echo_amp,'k');
legend('T2 decay', 'optimized','CFA 120');
xlabel('Echo Index');
ylabel('Signal');
title('Signal');   

subplot(122)
plot(refocus_alphas, '-');
hold on;
plot(x0_oneCH(2:end));
legend('optimized','x0');
ylim([0,190]);
title('RF pulse');
xlabel('Echo Index');
ylabel('Flip angle (deg)');
