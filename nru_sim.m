% L2 Simulation of the NRU model. 
% Prob of transmission calculated the ratio (total num of transmissions)/(total num of transmissions + total num of non-transmissions)


clear all;
close all;

W0 = 16;                            %backoff window size
m = 6;                              %from TS 37.213                        
K = 1;                              %retry limit
numUE = 1;                          %number of SL UEs
max_tt = 100000;                        %total number of transmission times

W = W0;
Wmax = 2^(m)*W0;                                %max backoff window size

succUE = 0;                                     %UE number that successfully sent packet
collUE = zeros(1, numUE);                       %collisions (set to 1) in previous iteration
W_UE = ones(1, numUE)*W0;                       %initial window size
stageUE = zeros(1, numUE);                      %re-transmission stage

cnt_notrans = zeros(1, numUE);                  %none-transmitting SL slots
cnt_trans = zeros(1, numUE);                    %transmissions by each UE


%initialization
bc = ones(1, numUE);
tt = 0;
transUE = [1:numUE];                            %indices of transmit UEs
notransUE = setdiff([1:numUE], transUE);        %indices of non-transmit UEs


while (tt < max_tt)

    tt = tt + 1;

    %A transmission occured in the previous iteration. Start a new backoff
    %process for the UEs that transmitted

    %transUE -- indices of previously transmitted UEs
    if(isempty(transUE) == 1)
        fprintf(1, 'Error in ind_transUE \n')
    end

    %% select a new backoff window/backoff value after a transmission

    if(length(transUE) == 1)
        %transmission is successful. Go to stage zero
        succUE = transUE;
        W_UE(succUE) = W0;
        stageUE(succUE) = 0;
        bc(succUE) = floor(W_UE(succUE)*rand);  %uniform in [0, W0-1]

    else
        % a collision occurred
        for uu = transUE

            if(stageUE(uu) < (m+K))
                stageUE(uu) = stageUE(uu) + 1;
                W_UE(uu) = min(Wmax, (2^stageUE(uu))*W0);
            else
                stageUE(uu) = 0;
                W_UE(uu) = W0;
            end

            bc(uu) = floor(W_UE(uu)*rand);  %uniform in [0, W_UE-1]

        end
    end


    %% update the backoff counters
    slot2Trans = min(bc);
    transUE = find(bc == slot2Trans);
    for uu = 1:numUE
        cnt_notrans(uu) = cnt_notrans(uu) + slot2Trans;
    end 
    for uu = 1:numUE

        slot_var = slot2Trans - bc(uu);
        if(slot_var == 0)

            %update the backoff counter
            bc(uu) = 0;
            cnt_trans(uu) = cnt_trans(uu) + 1;

        else
            %non-transmitting UEs
            bc_next = bc(uu) - slot2Trans;
            bcuu_temp = bc(uu);
            bc(uu) = bc_next;
        end


    end


end

% results

fprintf(1, '\nProbability of transmission for NR, numUE = %d  \n', numUE)
fprintf(1, 'tau_nr = %f \n', mean(cnt_trans) / ( mean(cnt_notrans)+mean(cnt_trans)));


