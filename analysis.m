% Analysis model for the SL-U only networks and NR-U only networks
clear all; close all;

%% Parameters from TS 37.213. Use p = 4 
W = 16;                         %TS 37.213. Table 4.2.1-1, p=4

%SL parameters
m_sl = 6;
e_sl = 1;                       %TS37.213 Clause 4.1.4.3 (K =  1,2, ..8)
num_sl = 1;
SL_slot = [250 500 1000];      %SL slot size (us)
% SL_slot = [500];

sensing_slot = 9;   %sensing slot size (us)
S = ceil(SL_slot/sensing_slot);
K = floor(6000./SL_slot); 
Tmcot = 6000;
r_SLU = [51.6, 49.7, 46.8];
for ii = 1:length(S)
    fun = @(x)tau_sl_rev_jan13(x, W, m_sl, e_sl, S(ii), num_sl);
    x0 = [0.4 0.4];
    options = optimoptions('fsolve', 'Display','none');
    [xvec, fval] = fsolve(fun, x0, options);
    if(norm(fval)> 1e-5)
        warning('Error in computation');
    end
    tau_sl(ii) = xvec(1);
    p_sl(ii) = xvec(2);


    %prob of successful transmission for SL-U network
    P_s_SL(ii) = num_sl * tau_sl(ii)*(1-tau_sl(ii))^(num_sl-1);

    %prob of collision for the SL-U network
    P_tr_SL(ii) = 1 - (1-tau_sl(ii))^num_sl;
    % P_c_SL(ii) = P_tr_SL(ii) - P_s_SL(ii);
    P_s_SL(ii) = 1/P_tr_SL(ii)*num_sl*tau_sl(ii)*(1-tau_sl(ii))^(num_sl-1);
    
%     Ts_SL = K(ii)* sensing_slot;
    Ts_SL = Tmcot*13/14;
    Tc_SL = Ts_SL;
    r_SL = r_SLU(ii);
    T_E = (1 - P_tr_SL(ii))*sensing_slot/S(ii) + P_tr_SL(ii)*P_s_SL(ii)*Ts_SL/S(ii) + P_tr_SL(ii)*(1- P_s_SL(ii))*Tc_SL/S(ii) + (S(ii) - 1)*sensing_slot/S(ii);

    Thpt = P_tr_SL(ii)*P_s_SL(ii)*Ts_SL*r_SL/T_E/S(ii);
    SL_Thpt(ii) = Thpt;
    
    fprintf(1, 'slot size (us) = %5d, number of nodes = %5d\n', ...
        SL_slot(ii), num_sl)
    fprintf(1, '\tProbability of transmission (tau_sl) = %f \n', tau_sl(ii))

    fprintf(1, '\tThroughput for SL = %f\n\n', Thpt);   
    
x = xvec;
[F, out] = tau_sl_rev_jan13(x, W, m_sl, e_sl, S(ii), num_sl);
fprintf(1, '\t S1 = %f, S2 = %f, S3 = %f \n', ...
    out.S1, out.S2, out.S3);
    

end

%compute the NR-U values
m_nr = m_sl;
e_nr = e_sl;
num_nr = num_sl;
fun = @(x)nru_tau(x, W, m_nr, e_nr, num_nr);
x0 = [0.05 0.1];
options = optimoptions('fsolve', 'Display','none');
[xvec, fval] = fsolve(fun, x0, options);
if(norm(fval)> 1e-5)
    warning('Error in computation');
end
tau_nr = xvec(1);
p_nr = xvec(2);


fprintf(1, 'tau_nr = %4.3f\n', tau_nr)
fprintf(1, 'Ratio tau_sl/tau_nr, tau_sl \n')
for ii = 1:length(tau_sl)
    fprintf(1, ['\tslot size = %4d, tau_sl/tau_nr = %4.3f,' ...
        '  tau_sl = %4.3f\n'], SL_slot(ii), tau_sl(ii)/tau_nr, tau_sl(ii))
end


m_nr = m_sl;
e_nr = e_sl;
num_nr = num_sl;
fun = @(x)nru_tau(x, W, m_nr, e_nr, num_nr);
x0 = [0.05 0.1];
options = optimoptions('fsolve', 'Display','none');
[xvec, fval] = fsolve(fun, x0, options);
if(norm(fval)> 1e-5)
    warning('Error in computation');
end
tau_nr = xvec(1);
p_nr = xvec(2);
%prob of successful transmission for NR-U network
P_s_NR  = num_nr * tau_nr *(1-tau_nr )^(num_nr-1);

%prob of collision for the NR-U network
P_tr_NR  = 1 - (1-tau_nr )^num_nr;
P_c_NR  = P_tr_NR  - P_s_NR ;
P_s_NR  = P_s_NR /P_tr_NR ;

% Ts_NR = 79*13/14;                     %fixed by Hao
Ts_NR = Tmcot*13/14;
Tc_NR = Ts_NR;
r_NR = 54.9;
T_E = (1 - P_tr_NR )*sensing_slot + P_tr_NR *P_s_NR *Ts_NR + P_c_NR *Tc_NR;
Thpt = P_tr_NR *P_s_NR *Ts_NR*r_NR/T_E;
NR_Thpt = Thpt;
fprintf(1, '\tThroughput for NR = %f\n\n', Thpt);   

%% NR-U equations (without coexistence)
function F =  nru_tau(x, W, m_nr, e_nr, num_nr)

tau_nr = x(1);
p_nr = x(2);

%NR system
var1 = (1-(2*p_nr)^(m_nr+1))/(1-2*p_nr);
var2 = (2^m_nr)*(p_nr^(m_nr+1) - p_nr^(m_nr+e_nr+1))/(1-p_nr);
var3 = (1-p_nr^(m_nr+e_nr+1))/2/(1-p_nr);
D1 = W/2*(var1+var2)+var3;

b00 = 1/(D1);

tau_n = (1-p_nr^(m_nr+e_nr+1))/(1-p_nr)*b00;
p_n = 1 - (1-tau_nr)^(num_nr-1);

F(1) = x(1) - tau_n;
F(2) = x(2) - p_n;

end


%% Calculate new tau_sl value

function [F, out] = tau_sl_rev_jan13(x, W, m, K, S, num_sl)

tau_sl = x(1);
p_sl = x(2); 
pb = p_sl;


% S1 value
var1 = (1-(2*p_sl)^(m+1))/(1-2*p_sl);
var2 = (2^m)*(p_sl^(m+1) - p_sl^(m+K+1))/(1-p_sl);
var3 = (1-p_sl^(m+K+1))/2/(1-p_sl);
S1 = W/2*(var1+var2)+var3;

% S2 and S3 values
var4 = (1 -(p_sl/2)^(m+1))/(1 - p_sl/2);
var5 = (p_sl^(m+1) - p_sl^(m+K+1))/(2^m*(1-p_sl));
S2 = var4/W + var5/W;
S2 = 1/W*(-var4 -var5) + 2*var3;

S3 = (1 - p_sl^(m+K+1))/(1-p_sl);


% F calculation
Fval = ((S-1)/2*((S+2*(1-pb)))+1)/(S*(1-pb/S));

% % b00 calculation
b00 = 1/(S1 + (Fval-1)*S2);
% 
% % tau calculation
% tau1 = (1-p_sl^(m+K+1))/(1-p_sl)*b00;


% % b00 calculation
% b00 = 1/(S1 + (S-1)/2*S3 - (S-1)*S2);


% tau calculation
tau1 = S3*b00;


%factor to account for the average number of sensing slots in an SL slot
p_t = (1 - (1 - tau_sl)^num_sl);

fac1 = S + S*p_t*log(p_t)/(1 - p_t)/2; 

tau_s = tau1*fac1;

p_s = 1 - (1-tau_sl)^(num_sl-1);

F(1) = x(1) - tau_s;
F(2) = x(2) - p_s;

out.S1 = S1;
out.S2 = S2;
out.S3 = S3;
out.b00 = b00;
out.tau1 = tau1;
% out.tau0 = tau0;
out.fac1 = fac1;
out.ps = p_s;
out.Fval = Fval;

end











