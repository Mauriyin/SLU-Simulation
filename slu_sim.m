% Simulation of the SLU model. 
% Time-continuous effect considered using tot_SLslot. Prob of transmission
% calculated the ratio (total num of transmissions)/(total num of SLslots)


clear all;
close all;

W0 = 16;                            %backoff window size
m = 6;                              %from TS 37.213                        
K = 1;                              %retry limit
SL_slottime = 1000;                  %SL-U slot size in mu secs [250, 500, 1000]
sensing_slottime = 9;               %sensing slot size in mu secs
numUE = 8;                          %number of SL UEs

S = ceil(SL_slottime/sensing_slottime); %sensing slots in a SL time slot (S value)
max_tt = 100000;                        %total number of transmission times

W = W0;
Wmax = 2^(m)*W0;                                %max backoff window size

succUE = 0;                                     %UE number that successfully sent packet
collUE = zeros(1, numUE);                       %collisions (set to 1) in previous iteration
W_UE = ones(1, numUE)*W0;                       %initial window size
stageUE = zeros(1, numUE);                      %re-transmission stage

cnt_bc = zeros(numUE, Wmax);                    %count backoff values
cnt_bs = zeros(numUE, S);                       %count gap states
cnt_bs0 = zeros(numUE, S);                      %count gap states when bc = 0
cnt_bs_S = zeros(1, numUE);                     %gap state 'S' transitioning from backoff state 2
cnt_notrans = zeros(1, numUE);                  %none-transmitting SL slots
cnt_SLslot = zeros(1, numUE);                   %total number of SL slots
cnt_trans = zeros(1, numUE);                    %transmissions by each UE
cnt_trans_other = zeros(1, numUE);              %other UEs transmit but not particular UE


%initialization
bc = ones(1, numUE);
tt = 0;
transUE = [1:numUE];                            %indices of transmit UEs
notransUE = setdiff([1:numUE], transUE);        %indices of non-transmit UEs
end_trans = 0;


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
        cnt_bc(succUE,bc(succUE)+1) = cnt_bc(succUE,bc(succUE)+1) + 1;

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
            cnt_bc(uu,bc(uu)+1) = cnt_bc(uu,bc(uu)+1) + 1;

        end
    end

    %UE uu is not transmitting but other UEs are transmitting
    for uu = 1:numUE
        ind2 = find(uu == transUE);
        if(isempty(ind2) == 1)
            cnt_trans_other(uu) = cnt_trans_other(uu) + 1;
        end
    end

    %% First SL slot boundary where a transmission occurs

    start_trans = end_trans;                            %range is [0, S-1] 
    [bc_min, ind_min] = min(bc);
    time_bcmin_zero = start_trans + bc_min;             %sensing slot when bc_min = 0

    SLslot = ceil(time_bcmin_zero/S);                   %SL slot boundary where a transmission occurs
    gap_len = SLslot*S - (start_trans + bc);
    slots2SLslot = SLslot*S - start_trans;              %sensing slots to transmission SL boundary slot

    notransUE_prev = notransUE;                         %none-transmitting UEs in the prev iteration
    transUE = find(gap_len >= 0);                       %transmitting UEs in the first SL slot
    notransUE = find(gap_len <0);                      
    %non-transmitting UEs in this SL slot while there is a transmission in trans_SLslot

    %% count no transmissions at SL slot boundaries

    if(start_trans == 0)
        missSLslot_allUE = SLslot;
        % SL boundary slots: [0, 1, 2, ..., SLslot]. First transmission
        % occurs at the SL slot SLslot

        %count the zero slot as the first missed SL bouudary slot only in the
        %UEs that transmitted immediately before this time
        ind1 = setdiff([1:numUE], notransUE_prev);          %UEs that transmitted immediately before
        cnt_notrans(ind1) = cnt_notrans(ind1) + missSLslot_allUE;
        cnt_notrans(notransUE_prev) = cnt_notrans(notransUE_prev) + missSLslot_allUE - 1;
        %UEs that did not transmit immediately before, cannot transmit at
        %this time

        %these UEs do not transmit at the SL transmission slot
        cnt_notrans(notransUE) = cnt_notrans(notransUE) + 1;

        %if a transmission occured in previous slot then count the slot
        %where transmission ended
        cnt_SLslot(ind1) = cnt_SLslot(ind1) + SLslot + 1;
        %if there were no transmissions in the previous slot then don't
        %count the slot where transmission ended
        cnt_SLslot(notransUE_prev) = cnt_SLslot(notransUE_prev) + SLslot;


    else
        missSLslot_allUE = SLslot - 1;
        %sensing slot zero, which is a SL boundary slot, is excluded

        cnt_notrans(notransUE) = cnt_notrans(notransUE) + missSLslot_allUE+1;
        cnt_notrans(transUE) = cnt_notrans(transUE) + missSLslot_allUE;

%         ind1 = setdiff([1:numUE], notransUE_prev);          %UEs that transmitted immediately before

        %if a transmission occured in previous slot then count the slot
        %where transmission ended
%         cnt_SLslot(ind1) = cnt_SLslot(ind1) + SLslot;
        %if there were no transmissions in the previous slot then don't
        %count the slot where transmission ended
%         cnt_SLslot(notransUE_prev) = cnt_SLslot(notransUE_prev) + SLslot;
        cnt_SLslot(1:numUE) = cnt_SLslot(1:numUE) + SLslot;

    end


    %% update the backoff counters at SL slot boundary

    for uu = 1:numUE

        slot_var = slots2SLslot - bc(uu);
        if(slot_var >= 0)
            % transmission occurs

            %update counting backoff values [0, bc(uu)-1]
            %if the UE transmitted in the previous slot, bc(uu) is included
            %in cnt_bc earlier. if the UE did not transmit in the previous
            %slot, bc(uu) is not included in this counting
            cnt_bc(uu, 1:bc(uu)) = cnt_bc(uu, 1:bc(uu)) + 1;
            

            % update gap states
            if(bc(uu)>= 1)
                gap_max = slots2SLslot - (bc(uu));
                %maximum size of gap occurs when bc(uu) reaches one
                cnt_bs(uu,1:gap_max) = cnt_bs(uu,1:gap_max) + 1;

            elseif(bc(uu) == 0)
                %if bc(uu) is zero send it via the same link as the bc = 1
                %path
                gap_max = slots2SLslot;
                cnt_bs0(uu,1:gap_max) = cnt_bs0(uu,1:gap_max) + 1;
            end

            %update the backoff counter
            bc(uu) = 0;
            cnt_trans(uu) = cnt_trans(uu) + 1;

        else
            %non-transmitting UEs
            bc_next = bc(uu) - slots2SLslot;
            cnt_bc(uu, bc_next+1:bc(uu)) = cnt_bc(uu, bc_next+1:bc(uu)) + 1;
            bcuu_temp = bc(uu);
            bc(uu) = bc_next;

            if(bc(uu) == 1)
                cnt_bs_S(uu) = cnt_bs_S(uu) + 1;
            end

        end


    end

    %select a random ending location 
    end_trans = floor(S*rand);                  %random end point in [0, S-1] 

end

%% results

fprintf(1, '\nProbability of transmission for S = %d, numUE = %d  \n', S, numUE)

%number of transmissions 
tau1 = mean(cnt_bc(:,1));

fprintf(1, '\tAverage transmissions:')
fprintf(1, 'Simulations = %7.2f, \t analysis = %7.2f \n \n', mean(cnt_trans), tau1);

tau0 = mean(cnt_notrans);
pr_trans_simu = tau1./(tau1+tau0);
t0 = zeros(1, numUE);
for uu = 1:numUE
    t0(uu) = sum(cnt_bc(uu,3:end));
end
pr_trans_cntstate = tau1/(tau1 + mean(t0)/S + mean(cnt_bs_S) + mean(cnt_bs(:,S)));
pr_trans_cntstate0 = tau1/(tau1 + mean(t0)/S + mean(cnt_bs_S));

fprintf(1, '\tSimultaions = %f, Counting states = %f, Cnt states0 = %f \n \n', ...
    pr_trans_simu, pr_trans_cntstate, pr_trans_cntstate0);


pr_trans2cntSLslot = tau1/(mean(cnt_SLslot));
fprintf(1, '\t pr(trans2cntSLslot) = %f \n', pr_trans2cntSLslot);


