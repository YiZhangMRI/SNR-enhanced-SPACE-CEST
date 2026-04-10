function [echo_amp,F_over_time,FZ_over_time] = func_forward_EPG(T1, T2, tesp, alpha_array, phi_array)
% forward extended phase graph (EPG) algorithm
% 
% Input
%   T1  
%   T2
%   tesp        - echo spacing, T1,T2,tesp should use same unit
%   alpha_array - flip angle array, unit is degree
%   phi_array   - phase of flip angle
% 
% Output
%   ehco_amp
%   F_over_time - states over time, 2nd dim is time
% 
% Example
%     T1 = 0.14;
%     T2 = 0.09;
%     tesp = 3.6/1000;
%     alpha_array = 120*ones(1,100);
%     echo_amp = func_forward_EPG(T1, T2, tesp, alpha_array);
%     figure;
%     plot(abs(echo_amp), 'r-');
% 
% Credits
%   https://github.com/matthias-weigel/EPG
%   https://github.com/mribri999/MRSignalsSeqs
% 
% See also bloch_simu_FSE

% TODO, use https://asciiflow.com/ to draw seq timing
% 
% Xingwang Yong, 20220108

%% input
num_states = 2*numel(alpha_array)+1;
F = zeros(3, num_states);
F(1,1) = 1;F(2,1) = 1; % init

if nargin < 5
    phi_array = zeros(size(alpha_array));
end


%% EPG simu
alpha_array = deg2rad(alpha_array);
E = R(T1, T2, tesp);

num_RF_pulse = numel(alpha_array);
echo_amp = zeros(1, num_RF_pulse);

all_z = zeros(1, num_RF_pulse); % for debug
FZ_sum = zeros(1, num_RF_pulse); 
FZ_square_sum = zeros(1, num_RF_pulse); 
F_over_time = zeros(numel(F)*2/3,num_RF_pulse);
FZ_over_time = zeros(numel(F)/3,num_RF_pulse);
for k = 1:num_RF_pulse       
%           |                       |                       |
%           |        1st iter       |       2nd iter        |
%           |<--------------------->|<--------------------->|
%           |                       |                       |
%           V                       V                       V
% 
%      (init state)              1st echo                 2nd echo
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           90         180		   echo		   180        echo         180         echo
%         				|			^			|                       |           	
%         			    |  		   / \			|		    ^           |
%         	|			|		  /   \			|		   / \          |           ^
%         	|			|		 /     \ 		|		  /	  \         |          / \
% ------------------------------------------------------------------------------------------------> time
%         	|			|						|			
%         	| S_op R_op | R_op S_op	  S_op R_op	|			
%         	|---------->|---------->|---------->|			
%         	|			|						|			
%                      T_op                    T_op
%               esp/2                          esp
%         	 <---------->			 <----------------------> 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = S_op(F); 
    F = R_op(E,F);
    F = T_op(F, alpha_array(k), phi_array(k));
    F = R_op(E,F);
    F = S_op(F);
    
    echo_amp(k) = F(1,1);
    
    all_z(k) = F(3,2);
    FZ_sum(k) = abs(sum(F, 'all'));
    FZ_square_sum(k) = abs(sum(F.^2, 'all'));
    
    
    F_over_time(:,k) = [F(1,end:-1:1),F(2,:)].'; % F(+k_max), F(+k_max-1), ..., F(0), F(0), F(-1), F(2), ..., F(-kmax)
    FZ_over_time(:,k) = F(3,end:-1:1).';
end

echo_amp = abs(echo_amp);
F_over_time(size(F_over_time,1)/2,:) = []; % delete the repeated F(0)
F_over_time = F_over_time(1:2:end,:); % only the primary pathway
FZ_over_time = FZ_over_time(2:2:end,:);% only the primary pathway, don't know why, copied from matthias-weigel
end


%% sub functions
function Fnew = S_op(F)
% Shift operator, shift over half of the echo spacing
% F is 3-by-N state matrix

Fnew = F;

Fnew(1,:) = circshift(F(1,:), 1);
Fnew(2,:) = circshift(F(2,:), -1);
Fnew(1,1) = conj(Fnew(2,1));
Fnew(2,end) = 0;
end

function Fnew = T_op(F, alpha, phi)
% RF transition operator
% F is 3-by-N state matrix

T = [        cos(alpha/2).^2        , exp(2*1i*phi)*sin(alpha/2).^2, -1i*exp(1i*phi)*sin(alpha);
     exp(-2*1i*phi)*sin(alpha/2).^2 ,        cos(alpha/2).^2       , 1i*exp(-1i*phi)*sin(alpha);
     -1i*0.5*exp(-1i*phi)*sin(alpha), 1i*0.5*exp(1i*phi)*sin(alpha),         cos(alpha)        ];
 
Fnew = T*F;
end

function Fnew = R_op(E, F)
% Relaxation operator
% F is 3-by-N state matrix

Fnew = E*F;
E1 = E(3,3);
Fnew(3,1) = Fnew(3,1) + 1-E1;
end

function E = R(T1, T2, tesp)
% Generate relaxation matrix

if T1==0
    E1 = 1;
else
    E1 = exp(-tesp/2/T1);
end

if T2==0
    E2 = 1;
else
    E2 = exp(-tesp/2/T2);
end

E = diag([E2,E2,E1]);
end