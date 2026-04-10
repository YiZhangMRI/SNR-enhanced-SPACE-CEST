function [obj, grad, FF] = SNR_res_obj(alpha_array, T1,T2,ESP,W,SNR_obj_weight,res_obj_weight,B1map,klim)
% copied from obj_EPG13.m, Alessandro Sbrizzi et al, MRM 2016


% ======================Added by Xingwang Yong======================
B1 = B1map;
[~, Nch] = size(B1);
numpulse = numel(alpha_array)/2/Nch;
params = alpha_array;  % real, imag
% params = repmat(alpha_array, 1,8);  % real, imag
frequencies = ones(1, numpulse);
c = W;
% if (numel(c)-2)~= (numel(alpha_array)/2-1)
%     error('Size mismatch');
% end
Ns = size(B1,1);
if Ns < 2
    error('Number of spatial points should be at least 2');
end
tmp_target = zeros(1, numpulse+1);
target = repmat(tmp_target(:),[ 1 Ns]); % target for all spatial locations
% ======================Added by Xingwang Yong======================





[nt] = size(target,1); 
[Ns, Nch] = size(B1);
np = length(params)/(2*Nch);
x = reshape(params(1:np*Nch),np,Nch).'; %real parts
y = reshape(params(np*Nch+1:end),np,Nch).'; % imaginary parts
params0 = x +1i*y;
params1 = zeros(Nch,sum(frequencies));
sumfreq = cumsum(frequencies);
for j = (length(frequencies)):-1:2
    params1(:,sumfreq(j-1)+1:sumfreq(j)) = repmat(params0(:,j),[1 frequencies(j)]);
end
params1(:,1:frequencies(1)) = repmat(params0(:,1),[1 frequencies(1)]);
params1 = params1.';
params2 = [real(params1(:).') , imag(params1(:).')];
params = params2;
%%
[Ns, Nch] = size(B1);
F0 = 0;

np = length(params)/(2*Nch);

x = reshape(params(1:np*Nch),np,Nch).'; %real parts
y = reshape(params(np*Nch+1:end),np,Nch).'; % imaginary parts
th = zeros(Ns,np,Nch);
for ch = 1:Nch
    th(:,:,ch) = B1(:,ch)*(x(ch,:)+1i*y(ch,:)); % effective theta
end
th = sum(th,3);

%th(:,2:end) = th(:,2:end)*exp(1i*pi/2); % not working!!
alph = abs(th);% effective alpha
ph = angle(th); % effective phi

% add CPMG phase
ph(:,2:end) = ph(:,2:end) + pi/2; % working!
xeff = real(th); % effective real
yeff = imag(th); % effective imag


kmax = 2*np - 1; % up to just before the next RF pulse
kmax = min([2*np-1, 2*klim-1]);
N = 3*(kmax-1)/2; % number of states in total
% split magnitude and phase

Npathway = inf;
klimit=false;
flipback=false;

% enforce pathway limit
if (N>Npathway)&&~klimit
    N=Npathway;
end

if klimit
    nr = size(th,2)-1; % number of refocus pulses
    kmax = nr - 1 + mod(nr,2);
    N = 1.5*(nr+mod(nr,2)); % number of states in total
    if mod(nr,2) 
        % odd
        KMAX = [1:2:kmax (kmax-2):-2:1];
    else
        %even
        KMAX = [1:2:kmax kmax:-2:1];
    end
    NMAX = 1.5*(KMAX+1);
else
    % previous implementation
    NMAX = 3:3:3*(np-1);
    NMAX(NMAX>N)=(N-mod(N,3));
end
    
    
%% ==== build Shift matrix, S with indices 
S = zeros([N N]);
%%% F(k>1) look @ states just BEFORE np+1 pulse
kidx = 4:3:N; % miss out F1+
sidx = kidx-3;
%idx = sub2ind([N N],kidx,sidx);
idx = kidx + N*(sidx-1);
S(idx)=1;

%%% F(k<1) 
kidx = 2:3:N;
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+3;
% ix = sub2ind([N N],kidx,sidx);
ix = kidx + N*(sidx-1);
S(ix)=1;

%%% Z states
kidx = 3:3:N;
%ix = sub2ind([N N],kidx,kidx);
ix = kidx + N*(kidx-1);
S(ix)=1;

%%% finally F1+ - relates to F-1-
S(1,2)=1;
%S(end-1,end-2) = 1; % not do this, it influences the truncated version!
S = sparse(S);
%% Relaxation =====: note, Z0 regrowth not handled (No Z0 state here)
E1=exp(-ESP/T1);
E2=exp(-ESP/T2);
%R = diag(repmat([E2 E2 E1],[1 kmax+1]));
R=eye(N);
ii = 1:(N/3);
R(3*N*(ii-1)+3*(ii-1)+1)=E2;
R(3*N*ii-2*N+3*(ii-1)+2)=E2;
R(3*N*ii-N+3*(ii-1)+3)=E1;

%%%% composites
RS=sparse(R*S);

%% F matrix (many elements zero, not efficient)
FF = zeros(N*2,np+1,Ns);

% Now the other states
for ns = 1:Ns % loop over space
    
    F = zeros([N*2 np+1]); %% records state 
    % initial state:
    F(3,1) = 1; 
    for jj=2:np+1 %loop over time

        if jj == 2 % excitation
            % Excitation pulse
            A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
            F([1:3 N+1:N+3],jj) = [real(A); imag(A)]*F(1:3,jj-1); %<---- state straight after excitation [F0 F0* Z0]
        elseif jj == 3 % first refocusing
            A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
            T = build_T_matrix_sub(A,3);
            F(1,jj) = exp(-0.5*ESP/T2)*F(1,jj-1);
            F(N+1,jj) = exp(-0.5*ESP/T2)*F(N+1,jj-1);
            F([1:3 N+1:N+3],jj) = [real(T), -imag(T);imag(T),real(T)]*[F(1:3,jj);F(N+1:N+3,jj)];  % Note by Xingwang Yong, this is not the signal, this is F_{k=-1}() states after refocus pulse, which will turn into signal later. the signal calculated this way would have a tesp/2 "displacement" with true signal. However, this is acceptable  
        else
            A = Trot_fun(alph(ns,jj-1),ph(ns,jj-1));
            temp1 = RS*F(1:N,jj-1);
            temp2 = RS*F(N+1:2*N,jj-1);
            temp2(1) = -temp2(1);% D*temp;
            F(:,jj) = [blockdiag_mult(real(A),temp1)-blockdiag_mult(imag(A),temp2); ...
                       blockdiag_mult(imag(A),temp1)+blockdiag_mult(real(A),temp2)];
        end
    end
    FF(:,:,ns) = F; 
end
% state = FF(2,:,:)+1i*FF(N+2,:,:); % not used, commented by Xingwang Yong
% state = diag(c)*(squeeze(state)-target); % not used, commented by Xingwang Yong
% obj = 0.5*norm(state(:))^2; % not used, commented by Xingwang Yong


% ======================Added by Xingwang Yong======================
state = FF(2,:,:)+1i*FF(N+2,:,:); 
state = diag(c)*(squeeze(state)-target); 
SNR_obj = -0.5*norm(state(:))^2; 

% res_obj = sum( (f(i)-f(i+1))^2 )
% sigs_diff = FF(2,3:end-1,:)+1i*FF(N+2,3:end-1,:) - FF(2,4:end,:) - 1i*FF(N+2,4:end,:);
% sigs_diff = squeeze(sigs_diff);
% res_obj = 0.5*norm(sigs_diff)^2;


sigs = squeeze(FF(2,3:end,:)+1i*FF(N+2,3:end,:));
s = 0;
% isOnlyOneSpatialPoint = isvector(sigs);
% if isOnlyOneSpatialPoint
%     sigs = sigs.';
% end
for iii=1:size(sigs,2) % loop over space
    sigs_one_pixel = sigs(:,iii);
    sigs_diff = repmat(sigs_one_pixel,1,numel(sigs_one_pixel)) - repmat(sigs_one_pixel.',numel(sigs_one_pixel),1);
%     s = s + norm(sigs_diff, 'fro')^2; % res_obj = sum( sum( (f(i)-f(j))^2 ) ), sum over all i,j

    sig_diff_norm = norm(sigs_diff,'fro')^2;
    sig_mean_norm = norm(sigs_one_pixel)^2;
    s = s + sig_diff_norm / sig_mean_norm; % res_obj_1 = sum( sum( (f(i)-f(j))^2 ) ), sum over all i,j; res_obj_2 = sum( fi*fi ); res_obj = res_obj_1/res_obj_2, this is somewhat similar to std()/mean()
end
res_obj = 0.5*s;


obj = SNR_obj_weight*SNR_obj + res_obj_weight*res_obj;
% ======================Added by Xingwang Yong======================



%% now adjoint states
if nargout > 1 % need gradient 

%C = zeros(2*N,2*N);
%C(2,2) = 1; % sampling matrix
%C(N+2,N+2) = 1; % sampling matrix

LL = zeros(np+1,N*2,Ns);
RRSS = [RS, RS*0;RS*0, RS];
RRSS(N+1,:) = -RRSS(N+1,:);
% Now the other states
for ns = 1:Ns % loop over space
    L = zeros(np+1,N*2);
    for jj=np:-1:1

        if jj == 1 % 
            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            T = build_T_matrix_sub(A,3);
            d = zeros(6,6);
            d(1) = exp(-0.5*ESP/T2);d(3+1,3+1) = exp(-0.5*ESP/T2);
            L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d+[0,c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns))),0,0,c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns))),0];
        
        elseif jj == np
        
            temp = zeros(1,2*N);
            temp(2) = c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns)));
            temp(N+2) = c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns)));
            L(jj,:) = temp; 
        else

            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            % NMAX is maximum state index
            kidx = 1:NMAX(end);
            temp = [blockdiag_mult(real(A).',L(jj+1,1:N).')+blockdiag_mult(imag(A).',L(jj+1,N+1:end).');...
                   -blockdiag_mult(imag(A).',L(jj+1,1:N).')+blockdiag_mult(real(A).',L(jj+1,N+1:end).')].';
            temp = temp*RRSS;
            temp1 = zeros(1,2*N);
            temp1(2) = c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns)));
            temp1(N+2) = c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns)));
            L(jj,:) = temp+temp1;%-(F(:,jj+1))'*C;

        end
    end
    LL(:,:,ns) = L; 
end

%% now gradient

GRAD = zeros(2*np*Nch,Ns);
for ns = 1:Ns % loop over space
    for jj = 2:np+1
        dT_daeff = dTrot_fun_da(alph(ns,jj-1),ph(ns,jj-1));
        dT_dpeff = dTrot_fun_dp(alph(ns,jj-1),ph(ns,jj-1));
        d = zeros(2*N,2*N);
        d(1,1) = exp(-0.5*ESP/T2);d(N+1,N+1) = exp(-0.5*ESP/T2);
        temp1 = d*FF(:,jj-1,ns);
        temp2 = RRSS*FF(:,jj-1,ns);
        TEMP1a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp1(1:N))-blockdiag_mult(imag(dT_daeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp1(1:N))+blockdiag_mult(real(dT_daeff),temp1(N+1:end))];
        TEMP1p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp1(1:N))-blockdiag_mult(imag(dT_dpeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp1(1:N))+blockdiag_mult(real(dT_dpeff),temp1(N+1:end))];
        TEMP2a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp2(1:N))-blockdiag_mult(imag(dT_daeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp2(1:N))+blockdiag_mult(real(dT_daeff),temp2(N+1:end))];
        TEMP2p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp2(1:N))-blockdiag_mult(imag(dT_dpeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp2(1:N))+blockdiag_mult(real(dT_dpeff),temp2(N+1:end))];
        
        for ch = 1:Nch
            dadxj = xeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch));
            dpdxj = xeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch))-yeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch));
            dadyj = -xeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch));
            dpdyj = xeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch));
            if jj == 3
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP1a+dpdxj*TEMP1p;
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP1a+dpdyj*TEMP1p;
            else
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP2a+dpdxj*TEMP2p; 
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP2a+dpdyj*TEMP2p; 
            end
        end
    end
end

grad = real(sum(GRAD,2));
grad = grad(:);
snr_obj_grad = -grad(:); % Added by Xingwang Yong, our SNR obj is obj=-1*sum(..), the original obj of optimal control EPG is obj=1*sum()


% % ==========Added by Xingwang Yong, calculate the adjoint state of resolution related obj, res_obj = sum( (f(i)-f(i+1))^2 )==========
% % (1) adjoint states
% LL = zeros(np+1,N*2,Ns);
% RRSS = [RS, RS*0;RS*0, RS];
% RRSS(N+1,:) = -RRSS(N+1,:);
% % Now the other states
% for ns = 1:Ns % loop over space
%     L = zeros(np+1,N*2);
%     for jj=np:-1:1
% 
%         if jj == 1 % 
%             A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
%             T = build_T_matrix_sub(A,3);
%             d = zeros(6,6);
%             d(1) = exp(-0.5*ESP/T2);d(3+1,3+1) = exp(-0.5*ESP/T2);
% %             L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d+[0,c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns))),0,0,c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns))),0];
%             L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d + ...
%                 [0, 2*FF(2,jj+1,ns)-FF(2,jj+2,ns)-FF(2,jj,ns), 0, 0, 2*FF(N+2,jj+1,ns)-FF(N+2,jj+2,ns)-FF(N+2,jj,ns), 0];
%         
%         elseif jj == np
%         
%             temp = zeros(1,2*N);
%             temp(2)   = FF(2,jj+1,ns) - FF(2,jj,ns);
%             temp(N+2) = FF(N+2,jj+1,ns) - FF(N+2,jj,ns);
%             L(jj,:) = temp; 
%         else
% 
%             A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
%             % NMAX is maximum state index
%             kidx = 1:NMAX(end);
%             temp = [blockdiag_mult(real(A).',L(jj+1,1:N).')+blockdiag_mult(imag(A).',L(jj+1,N+1:end).');...
%                    -blockdiag_mult(imag(A).',L(jj+1,1:N).')+blockdiag_mult(real(A).',L(jj+1,N+1:end).')].';
%             temp = temp*RRSS;
%             temp1 = zeros(1,2*N);                     
%             temp1(2)   = 2*FF(2,jj+1,ns) - FF(2,jj+2,ns) - FF(2,jj,ns);
%             temp1(N+2) = 2*FF(N+2,jj+1,ns) - FF(N+2,jj+2,ns) - FF(N+2,jj,ns);                                    
%             L(jj,:) = temp+temp1;%-(F(:,jj+1))'*C;            
% 
%         end
%     end
%     LL(:,:,ns) = L; 
% end
% 
% 
% 
% % (2) gradient
% GRAD = zeros(2*np*Nch,Ns);
% for ns = 1:Ns % loop over space
%     for jj = 2:np+1
%         dT_daeff = dTrot_fun_da(alph(ns,jj-1),ph(ns,jj-1));
%         dT_dpeff = dTrot_fun_dp(alph(ns,jj-1),ph(ns,jj-1));
%         d = zeros(2*N,2*N);
%         d(1,1) = exp(-0.5*ESP/T2);d(N+1,N+1) = exp(-0.5*ESP/T2);
%         temp1 = d*FF(:,jj-1,ns);
%         temp2 = RRSS*FF(:,jj-1,ns);
%         TEMP1a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp1(1:N))-blockdiag_mult(imag(dT_daeff),temp1(N+1:end));...
%                                 blockdiag_mult(imag(dT_daeff),temp1(1:N))+blockdiag_mult(real(dT_daeff),temp1(N+1:end))];
%         TEMP1p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp1(1:N))-blockdiag_mult(imag(dT_dpeff),temp1(N+1:end));...
%                                 blockdiag_mult(imag(dT_dpeff),temp1(1:N))+blockdiag_mult(real(dT_dpeff),temp1(N+1:end))];
%         TEMP2a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp2(1:N))-blockdiag_mult(imag(dT_daeff),temp2(N+1:end));...
%                                 blockdiag_mult(imag(dT_daeff),temp2(1:N))+blockdiag_mult(real(dT_daeff),temp2(N+1:end))];
%         TEMP2p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp2(1:N))-blockdiag_mult(imag(dT_dpeff),temp2(N+1:end));...
%                                 blockdiag_mult(imag(dT_dpeff),temp2(1:N))+blockdiag_mult(real(dT_dpeff),temp2(N+1:end))];
%         
%         for ch = 1:Nch
%             dadxj = xeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch));
%             dpdxj = xeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch))-yeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch));
%             dadyj = -xeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch));
%             dpdyj = xeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch));
%             if jj == 3
%                 GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP1a+dpdxj*TEMP1p;
%                 GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP1a+dpdyj*TEMP1p;
%             else
%                 GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP2a+dpdxj*TEMP2p; 
%                 GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP2a+dpdyj*TEMP2p; 
%             end
%         end
%     end
% end
% 
% grad2 = real(sum(GRAD,2));
% grad2 = grad2(:);
% 
% 
% 
% 
% grad = SNR_obj_weight*grad + res_obj_weight*grad2;  % the overall gradient
% % ==========Added by Xingwang Yong, calculate the adjoint state of resolution related obj, res_obj = sum( (f(i)-f(i+1))^2 )==========





% % ==========Added by Xingwang Yong, calculate the adjoint state of resolution related obj, res_obj = sum( sum( (f(i)-f(j))^2 ) ), sum over all i,j==========
% % (1) adjoint states
% LL = zeros(np+1,N*2,Ns);
% RRSS = [RS, RS*0;RS*0, RS];
% RRSS(N+1,:) = -RRSS(N+1,:);
% % Now the other states
% for ns = 1:Ns % loop over space
%     L = zeros(np+1,N*2);
%     for jj=np:-1:1
% 
%         if jj == 1 % 
%             A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
%             T = build_T_matrix_sub(A,3);
%             d = zeros(6,6);
%             d(1) = exp(-0.5*ESP/T2);d(3+1,3+1) = exp(-0.5*ESP/T2);
% %             L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d+[0,c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns))),0,0,c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns))),0];
%             L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d + ...
%                 [0, 2*sum( FF(2,jj+1,ns)-FF(2,3:end,ns) ), 0, 0, 2*sum(  FF(N+2,jj+1,ns)-FF(N+2,3:end,ns) ), 0]; 
%         
%         elseif jj == np
%         
%             temp = zeros(1,2*N);
%             temp(2)   = 2*sum(  FF(2,jj+1,ns)   - FF(2,3:end,ns)  );
%             temp(N+2) = 2*sum(  FF(N+2,jj+1,ns) - FF(N+2,3:end,ns)  );
%             L(jj,:) = temp; 
%         else
% 
%             A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
%             % NMAX is maximum state index
%             kidx = 1:NMAX(end);
%             temp = [blockdiag_mult(real(A).',L(jj+1,1:N).')+blockdiag_mult(imag(A).',L(jj+1,N+1:end).');...
%                    -blockdiag_mult(imag(A).',L(jj+1,1:N).')+blockdiag_mult(real(A).',L(jj+1,N+1:end).')].';
%             temp = temp*RRSS;
%             temp1 = zeros(1,2*N);  
%             temp1(2)   = 2*sum(  FF(2,jj+1,ns)   - FF(2,3:end,ns)  );
%             temp1(N+2) = 2*sum(  FF(N+2,jj+1,ns) - FF(N+2,3:end,ns)  );                                  
%             L(jj,:) = temp+temp1;%-(F(:,jj+1))'*C;            
% 
%         end
%     end
%     LL(:,:,ns) = L; 
% end
% 
% 
% 
% % (2) gradient
% GRAD = zeros(2*np*Nch,Ns);
% for ns = 1:Ns % loop over space
%     for jj = 2:np+1
%         dT_daeff = dTrot_fun_da(alph(ns,jj-1),ph(ns,jj-1));
%         dT_dpeff = dTrot_fun_dp(alph(ns,jj-1),ph(ns,jj-1));
%         d = zeros(2*N,2*N);
%         d(1,1) = exp(-0.5*ESP/T2);d(N+1,N+1) = exp(-0.5*ESP/T2);
%         temp1 = d*FF(:,jj-1,ns);
%         temp2 = RRSS*FF(:,jj-1,ns);
%         TEMP1a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp1(1:N))-blockdiag_mult(imag(dT_daeff),temp1(N+1:end));...
%                                 blockdiag_mult(imag(dT_daeff),temp1(1:N))+blockdiag_mult(real(dT_daeff),temp1(N+1:end))];
%         TEMP1p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp1(1:N))-blockdiag_mult(imag(dT_dpeff),temp1(N+1:end));...
%                                 blockdiag_mult(imag(dT_dpeff),temp1(1:N))+blockdiag_mult(real(dT_dpeff),temp1(N+1:end))];
%         TEMP2a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp2(1:N))-blockdiag_mult(imag(dT_daeff),temp2(N+1:end));...
%                                 blockdiag_mult(imag(dT_daeff),temp2(1:N))+blockdiag_mult(real(dT_daeff),temp2(N+1:end))];
%         TEMP2p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp2(1:N))-blockdiag_mult(imag(dT_dpeff),temp2(N+1:end));...
%                                 blockdiag_mult(imag(dT_dpeff),temp2(1:N))+blockdiag_mult(real(dT_dpeff),temp2(N+1:end))];
%         
%         for ch = 1:Nch
%             dadxj = xeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch));
%             dpdxj = xeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch))-yeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch));
%             dadyj = -xeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch));
%             dpdyj = xeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch));
%             if jj == 3
%                 GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP1a+dpdxj*TEMP1p;
%                 GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP1a+dpdyj*TEMP1p;
%             else
%                 GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP2a+dpdxj*TEMP2p; 
%                 GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP2a+dpdyj*TEMP2p; 
%             end
%         end
%     end
% end
% 
% grad2 = real(sum(GRAD,2));
% grad2 = grad2(:);
% 
% 
% 
% 
% grad = SNR_obj_weight*grad + res_obj_weight*grad2;  % the overall gradient
% % ==========Added by Xingwang Yong, calculate the adjoint state of resolution related obj, res_obj = sum( sum( (f(i)-f(j))^2 ) ), sum over all i,j==========





% ==========Added by Xingwang Yong, calculate the adjoint state of resolution related obj, res_obj_1 = sum( sum( (f(i)-f(j))^2 ) ), sum over all i,j; res_obj_2 = sum( fi*fi ); res_obj = res_obj_1/res_obj_2, this is somewhat similar to std()/mean()==========
% we want to calculate the derivative of g(x)/h(x), first we calculate the 
% derivative of g(x) and h(x) separately, then combine them


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>the denominator, h(x)>>>>>>>>>>>>>>>>>>>>>>
c(3:end) = 1; % here c is simply a "sampling mask", not the weight used in SNR objective
LL = zeros(np+1,N*2,Ns);
RRSS = [RS, RS*0;RS*0, RS];
RRSS(N+1,:) = -RRSS(N+1,:);
% Now the other states
for ns = 1:Ns % loop over space
    L = zeros(np+1,N*2);
    for jj=np:-1:1

        if jj == 1 % 
            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            T = build_T_matrix_sub(A,3);
            d = zeros(6,6);
            d(1) = exp(-0.5*ESP/T2);d(3+1,3+1) = exp(-0.5*ESP/T2);
            L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d+[0,c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns))),0,0,c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns))),0];
        
        elseif jj == np
        
            temp = zeros(1,2*N);
            temp(2) = c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns)));
            temp(N+2) = c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns)));
            L(jj,:) = temp; 
        else

            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            % NMAX is maximum state index
            kidx = 1:NMAX(end);
            temp = [blockdiag_mult(real(A).',L(jj+1,1:N).')+blockdiag_mult(imag(A).',L(jj+1,N+1:end).');...
                   -blockdiag_mult(imag(A).',L(jj+1,1:N).')+blockdiag_mult(real(A).',L(jj+1,N+1:end).')].';
            temp = temp*RRSS;
            temp1 = zeros(1,2*N);
            temp1(2) = c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns)));
            temp1(N+2) = c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns)));
            L(jj,:) = temp+temp1;%-(F(:,jj+1))'*C;

        end
    end
    LL(:,:,ns) = L; 
end

% now gradient
GRAD = zeros(2*np*Nch,Ns);
for ns = 1:Ns % loop over space
    for jj = 2:np+1
        dT_daeff = dTrot_fun_da(alph(ns,jj-1),ph(ns,jj-1));
        dT_dpeff = dTrot_fun_dp(alph(ns,jj-1),ph(ns,jj-1));
        d = zeros(2*N,2*N);
        d(1,1) = exp(-0.5*ESP/T2);d(N+1,N+1) = exp(-0.5*ESP/T2);
        temp1 = d*FF(:,jj-1,ns);
        temp2 = RRSS*FF(:,jj-1,ns);
        TEMP1a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp1(1:N))-blockdiag_mult(imag(dT_daeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp1(1:N))+blockdiag_mult(real(dT_daeff),temp1(N+1:end))];
        TEMP1p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp1(1:N))-blockdiag_mult(imag(dT_dpeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp1(1:N))+blockdiag_mult(real(dT_dpeff),temp1(N+1:end))];
        TEMP2a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp2(1:N))-blockdiag_mult(imag(dT_daeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp2(1:N))+blockdiag_mult(real(dT_daeff),temp2(N+1:end))];
        TEMP2p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp2(1:N))-blockdiag_mult(imag(dT_dpeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp2(1:N))+blockdiag_mult(real(dT_dpeff),temp2(N+1:end))];
        
        for ch = 1:Nch
            dadxj = xeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch));
            dpdxj = xeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch))-yeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch));
            dadyj = -xeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch));
            dpdyj = xeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch));
            if jj == 3
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP1a+dpdxj*TEMP1p;
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP1a+dpdyj*TEMP1p;
            else
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP2a+dpdxj*TEMP2p; 
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP2a+dpdyj*TEMP2p; 
            end
        end
    end
end

grad = real(sum(GRAD,2));
grad_denominator = grad(:);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<the denominator, h(x)<<<<<<<<<<<<<<<<<<<<<<


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>the nominator, g(x)>>>>>>>>>>>>>>>>>>>>>>>>
LL = zeros(np+1,N*2,Ns);
RRSS = [RS, RS*0;RS*0, RS];
RRSS(N+1,:) = -RRSS(N+1,:);
% Now the other states
for ns = 1:Ns % loop over space
    L = zeros(np+1,N*2);
    for jj=np:-1:1

        if jj == 1 % 
            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            T = build_T_matrix_sub(A,3);
            d = zeros(6,6);
            d(1) = exp(-0.5*ESP/T2);d(3+1,3+1) = exp(-0.5*ESP/T2);
%             L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d+[0,c(jj+1)*(FF(2,jj+1,ns)-real(target(jj+1,ns))),0,0,c(jj+1)*(FF(N+2,jj+1,ns)-imag(target(jj+1,ns))),0];
            L(jj,[1:3 N+1:N+3]) = L(jj+1,[1:3 N+1:N+3])*[real(T) -imag(T);imag(T) real(T)]*d + ...
                [0, 2*sum( FF(2,jj+1,ns)-FF(2,3:end,ns) ), 0, 0, 2*sum(  FF(N+2,jj+1,ns)-FF(N+2,3:end,ns) ), 0]; 
        
        elseif jj == np
        
            temp = zeros(1,2*N);
            temp(2)   = 2*sum(  FF(2,jj+1,ns)   - FF(2,3:end,ns)  );
            temp(N+2) = 2*sum(  FF(N+2,jj+1,ns) - FF(N+2,3:end,ns)  );
            L(jj,:) = temp; 
        else

            A = Trot_fun(alph(ns,jj+1),ph(ns,jj+1));
            % NMAX is maximum state index
            kidx = 1:NMAX(end);
            temp = [blockdiag_mult(real(A).',L(jj+1,1:N).')+blockdiag_mult(imag(A).',L(jj+1,N+1:end).');...
                   -blockdiag_mult(imag(A).',L(jj+1,1:N).')+blockdiag_mult(real(A).',L(jj+1,N+1:end).')].';
            temp = temp*RRSS;
            temp1 = zeros(1,2*N);  
            temp1(2)   = 2*sum(  FF(2,jj+1,ns)   - FF(2,3:end,ns)  );
            temp1(N+2) = 2*sum(  FF(N+2,jj+1,ns) - FF(N+2,3:end,ns)  );                                  
            L(jj,:) = temp+temp1;%-(F(:,jj+1))'*C;            

        end
    end
    LL(:,:,ns) = L; 
end

% gradient
GRAD = zeros(2*np*Nch,Ns);
for ns = 1:Ns % loop over space
    for jj = 2:np+1
        dT_daeff = dTrot_fun_da(alph(ns,jj-1),ph(ns,jj-1));
        dT_dpeff = dTrot_fun_dp(alph(ns,jj-1),ph(ns,jj-1));
        d = zeros(2*N,2*N);
        d(1,1) = exp(-0.5*ESP/T2);d(N+1,N+1) = exp(-0.5*ESP/T2);
        temp1 = d*FF(:,jj-1,ns);
        temp2 = RRSS*FF(:,jj-1,ns);
        TEMP1a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp1(1:N))-blockdiag_mult(imag(dT_daeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp1(1:N))+blockdiag_mult(real(dT_daeff),temp1(N+1:end))];
        TEMP1p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp1(1:N))-blockdiag_mult(imag(dT_dpeff),temp1(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp1(1:N))+blockdiag_mult(real(dT_dpeff),temp1(N+1:end))];
        TEMP2a = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_daeff),temp2(1:N))-blockdiag_mult(imag(dT_daeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_daeff),temp2(1:N))+blockdiag_mult(real(dT_daeff),temp2(N+1:end))];
        TEMP2p = LL(jj-1,:,ns)*[blockdiag_mult(real(dT_dpeff),temp2(1:N))-blockdiag_mult(imag(dT_dpeff),temp2(N+1:end));...
                                blockdiag_mult(imag(dT_dpeff),temp2(1:N))+blockdiag_mult(real(dT_dpeff),temp2(N+1:end))];
        
        for ch = 1:Nch
            dadxj = xeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch));
            dpdxj = xeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch))-yeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch));
            dadyj = -xeff(ns,jj-1)/alph(ns,jj-1)*imag(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)*real(B1(ns,ch));
            dpdyj = xeff(ns,jj-1)/alph(ns,jj-1)^2*real(B1(ns,ch))+yeff(ns,jj-1)/alph(ns,jj-1)^2*imag(B1(ns,ch));
            if jj == 3
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP1a+dpdxj*TEMP1p;
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP1a+dpdyj*TEMP1p;
            else
                GRAD(jj-1+(ch-1)*np,ns) = dadxj*TEMP2a+dpdxj*TEMP2p; 
                GRAD(np*Nch+jj-1+(ch-1)*np,ns) = dadyj*TEMP2a+dpdyj*TEMP2p; 
            end
        end
    end
end

grad2 = real(sum(GRAD,2));
grad_nominator = grad2(:);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<the nominator, g(x)<<<<<<<<<<<<<<<<<<<<<<<<


res_obj_grad = ( grad_nominator*sig_mean_norm - grad_denominator*sig_diff_norm ) / sig_mean_norm^2;


grad = SNR_obj_weight*snr_obj_grad + res_obj_weight*res_obj_grad;  % the overall gradient
% ==========Added by Xingwang Yong, calculate the adjoint state of resolution related obj, res_obj_1 = sum( sum( (f(i)-f(j))^2 ) ), sum over all i,j; res_obj_2 = sum( fi*fi ); res_obj = res_obj_1/res_obj_2, this is somewhat similar to std()/mean()==========




%% with respect to parameters
A = zeros(length(frequencies),sum(frequencies));
for j = 2:length(frequencies)
    A(j,sumfreq(j-1)+1:sumfreq(j-1)+frequencies(j)) = 1;
end
A(1,1:frequencies(1)) = 1;
A = sparse(A);
A = kron(eye(2*Nch),A);
grad = A*grad;


end % of if nargout>1



end % of function()



% reduced version
function T =  build_T_matrix_sub(AA,nn)
T=zeros([nn nn]);
%          for ii = 1:(nn/3);
%             T(3*nn*(ii-1)+3*(ii-1)+1)=AA(1);
%             T(3*nn*(ii-1)+3*(ii-1)+2)=AA(2);
%             T(3*nn*(ii-1)+3*(ii-1)+3)=AA(3);
%             T(3*nn*ii-2*nn+3*(ii-1)+1)=AA(4);
%             T(3*nn*ii-2*nn+3*(ii-1)+2)=AA(5);
%             T(3*nn*ii-2*nn+3*(ii-1)+3)=AA(6);
%             T(3*nn*ii-nn+3*(ii-1)+1)=AA(7);
%             T(3*nn*ii-nn+3*(ii-1)+2)=AA(8);
%             T(3*nn*ii-nn+3*(ii-1)+3)=AA(9);
%          end
ind = 1:(nn/3);
T(3*nn*(ind-1)+3*(ind-1)+1)=AA(1);
T(3*nn*(ind-1)+3*(ind-1)+2)=AA(2);
T(3*nn*(ind-1)+3*(ind-1)+3)=AA(3);
T(3*nn*ind-2*nn+3*(ind-1)+1)=AA(4);
T(3*nn*ind-2*nn+3*(ind-1)+2)=AA(5);
T(3*nn*ind-2*nn+3*(ind-1)+3)=AA(6);
T(3*nn*ind-nn+3*(ind-1)+1)=AA(7);
T(3*nn*ind-nn+3*(ind-1)+2)=AA(8);
T(3*nn*ind-nn+3*(ind-1)+3)=AA(9);

T = sparse(T);
%        T = sparse(kron(eye(nn/3),AA));% slow
end

function w = blockdiag_mult(A,v)
%given A square matrix, it returns w =  kron(A,eye(N))*v efficiently
nA = size(A,1); %square
V = reshape(v,nA,length(v)/nA);
W = A*V;
w = W(:);
end


% Rotation matrix direct definition
function T = Trot_fun(a,p)

T = zeros([3 3]);
T(1) = cos(a/2).^2;
T(2) = exp(-2*1i*p)*(sin(a/2)).^2;
T(3) = -0.5*1i*exp(-1i*p)*sin(a);
T(4) = conj(T(2));
T(5) = T(1);
T(6) = 0.5*1i*exp(1i*p)*sin(a);
T(7) = -1i*exp(1i*p)*sin(a);
T(8) = 1i*exp(-1i*p)*sin(a);
T(9) = cos(a);
end

function T = dTrot_fun_da(a,p)
% derivative Rotation matrix w.r.t. a
T = zeros([3 3]);
T(1) = -sin(a)/2;%cos(a/2).^2;
T(2) = 0.5*exp(-2*1i*p)*sin(a);%exp(-2*1i*p)*(sin(a/2)).^2;
T(3) = -0.5*1i*exp(-1i*p)*cos(a);%-0.5*1i*exp(-1i*p)*sin(a);
T(4) = conj(T(2));%conj(T(2));
T(5) = T(1);%T(1);
T(6) = 0.5*1i*exp(1i*p)*cos(a);%0.5*1i*exp(1i*p)*sin(a);
T(7) = -1i*exp(1i*p)*cos(a);%-1i*exp(1i*p)*sin(a);
T(8) = 1i*exp(-1i*p)*cos(a);%1i*exp(-1i*p)*sin(a);
T(9) = -sin(a);%cos(a);
end
function T = dTrot_fun_dp(a,p)
% derivative Rotation matrix w.r.t. p
T = zeros([3 3]);
T(1) = 0;%cos(a/2).^2;
T(2) = -2*1i*exp(-2*1i*p)*(sin(a/2)).^2;
T(3) = -0.5*exp(-1i*p)*sin(a);
T(4) = conj(T(2));
T(5) = T(1);%T(1);
T(6) = -0.5*exp(1i*p)*sin(a);
T(7) = exp(1i*p)*sin(a);
T(8) = exp(-1i*p)*sin(a);
T(9) = 0;%cos(a);
end
