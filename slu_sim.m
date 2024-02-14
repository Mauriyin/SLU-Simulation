% Description:  
%               Packet-level simulation of the SL-U model. SL-U network operating without
%               other networks. 802.11 DCF modified to account for the gap
%               duration and randomization of packet lengths
% 
% Input:       
%               W0 - Initial backoff window size
%               m  - Retransmission stage
%               K  - Retry limit. Try the maximum window size 2^m*(W0) for
%               K times
%               SL_slottime         - Sidelink slot size in us 
%               sensing_slottime    - sensing slot time in us
%               numUE               - Number of UEs in the SL-U network
%               max_tt              - Number of packet transmissions
%
% Output:
%               Transmission probability (tau_{SL})
%

clear all; 
close all;

% Input parameters
W0 = 16;                            %Initial backoff window size
m = 6;                              %From TS 37.213. Retransmission stage                        
K = 1;                              %From TS 37.213. Retry limit
SL_slottime = 1000;                 %SL-U slot size in mu secs [250, 500, 1000]
sensing_slottime = 9;               %Sensing slot size in mu secs (TS 37.213)
numUE = 8;                          %Number of SL UEs
max_tt = 100000;                    %Total number of transmission times (for all UEs)

% Derived values
S = ceil(SL_slottime/sensing_slottime);         %Sensing slots in a SL time slot
Wmax = 2^(m)*W0;                                %Max backoff window size


%%  Initialization
W_UE = ones(1, numUE)*W0;                       %initial window size
stageUE = zeros(1, numUE);                      %re-transmission stage
bc = ones(1, numUE);                            %backoff value
tt = 0;                                         %count packet transmission
transUE = [1:numUE];                            %indices of transmit UEs
notransUE = setdiff([1:numUE], transUE);        %indices of non-transmit UEs
end_trans = 0;                                  %radom position U[0, S-1] where each 
%                                               packet transmission ends

% Counting variables
cnt_bc = zeros(numUE, Wmax);                    %count backoff values
cnt_bs = zeros(numUE, S);                       %count gap state values
cnt_bs0 = zeros(numUE, S);                      %count gap states when bc = 0
cnt_bs_S = zeros(1, numUE);                     %gap state 'S' transitioning from backoff state 2
cnt_notrans = zeros(1, numUE);                  %none-transmitting SL slots
cnt_SLslot = zeros(1, numUE);                   %total number of SL slots

%% Iterate over each packet transmission

while (tt < max_tt)

    tt = tt + 1;                    %count number of packet transmissions

    %A transmission occured in the previous iteration. Start a new backoff
    %process for the UEs that transmitted

    %transUE -- indices of previously transmitted UEs
    if(isempty(transUE) == 1)
        fprintf(1, 'Error in transUE \n')
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

    %% First SL slot boundary where a transmission occurs

    start_trans = end_trans;                            %start sensing slot position U[0, S-1] 
    [bc_min, ind_min] = min(bc);
    time_bcmin_zero = start_trans + bc_min;             %sensing slot when bc_min = 0

    SLslot = ceil(time_bcmin_zero/S);                   %First SL slot boundary where a transmission occurs
    gap_len = SLslot*S - (start_trans + bc);
    slots2SLslot = SLslot*S - start_trans;              %Number of sensing slots to the  SL boundary slot

    notransUE_prev = notransUE;                         %none-transmitting UEs in the prev iteration
    transUE = find(gap_len >= 0);                       %transmitting UEs at the first SL slot
    notransUE = find(gap_len <0);                      
    %non-transmitting UEs in this SL slot while there is a transmission in SLslot

    %% count number of transmissions at SL slot boundaries

    if(start_trans == 0)
        missSLslot_allUE = SLslot;
        % SL boundary slots: [0, 1, 2, ..., SLslot]. First transmission
        % occurs at the SL slot SLslot

        %count the zero slot as the first missed SL boundary slot only in the
        %UEs that transmitted immediately before this time
        ind1 = setdiff([1:numUE], notransUE_prev);          %UEs that transmitted immediately before
        cnt_notrans(ind1) = cnt_notrans(ind1) + missSLslot_allUE;
        if(missSLslot_allUE ~= 0)
            cnt_notrans(notransUE_prev) = cnt_notrans(notransUE_prev) + missSLslot_allUE - 1;
            %UEs that did not transmit immediately before, cannot transmit at
            %this time
        end
 
        %these UEs do not transmit at the current SL transmission slot
        cnt_notrans(notransUE) = cnt_notrans(notransUE) + 1;

        %if a transmission occured in previous slot then count the slot
        %where that transmission ended. Slot # zero 
        cnt_SLslot(ind1) = cnt_SLslot(ind1) + SLslot + 1;
        %if there were no transmissions in the previous slot then don't
        %count the slot where transmission ended. (Count only the slots
        %where the UE could potentially transmit)
        cnt_SLslot(notransUE_prev) = cnt_SLslot(notransUE_prev) + SLslot;


    else
        missSLslot_allUE = SLslot - 1;
        %sensing slot zero, which is a SL boundary slot, is excluded
        % Effective SL slots are in range [1 SLslot]
        if(SLslot == 0)
            fprintf(1, 'Error in SLslot calculation \n')
        end

        cnt_notrans(notransUE) = cnt_notrans(notransUE) + missSLslot_allUE + 1;
        cnt_notrans(transUE) = cnt_notrans(transUE) + missSLslot_allUE;

        % Transmissions cannot take place in SLslot = 0 slot
        % Effective SL slots are in range [1 SLslot]
        cnt_SLslot(1:numUE) = cnt_SLslot(1:numUE) + SLslot;

    end


    %% update the backoff counters at SL slot boundary

    for uu = 1:numUE

        slot_var = slots2SLslot - bc(uu);
        if(slot_var >= 0)
            % Transmitting UEs

            %update counting backoff values [0, bc(uu)-1]
            cnt_bc(uu, 1:bc(uu)) = cnt_bc(uu, 1:bc(uu)) + 1;

            if(bc(uu)>= 1)
                %update gap state 1 path
                
                gap_max = slots2SLslot - (bc(uu));
                %maximum size of gap occurs when bc(uu) reaches one
                cnt_bs(uu,1:gap_max) = cnt_bs(uu,1:gap_max) + 1;

            elseif(bc(uu) == 0)
                %if bc(uu) is zero send it to the gap state 0 path
                gap_max = slots2SLslot;
                cnt_bs0(uu,1:gap_max) = cnt_bs0(uu,1:gap_max) + 1;
            end

            %update the backoff counter
            bc(uu) = 0;

        else
            %non-transmitting UEs
            bc_next = bc(uu) - slots2SLslot;
            cnt_bc(uu, bc_next+1:bc(uu)) = cnt_bc(uu, bc_next+1:bc(uu)) + 1;
            bcuu_temp = bc(uu);
            bc(uu) = bc_next;

            if(bc(uu) == 1)
                %gap state equal to S
                cnt_bs_S(uu) = cnt_bs_S(uu) + 1;
            end

        end


    end

    %select a random ending location for all packets
    end_trans = floor(S*rand);                  %random end point in [0, S-1] 

end

%% Results

tau1 = mean(cnt_bc(:,1));                           %number of transmissions 
prob_transmission = tau1/mean(cnt_SLslot);          %probability of transmission -- tau_{SL}


fprintf(1, '\nResults: SL-U slot size (us) = %5d,  S = %d, numUE = %d  \n', ...
    SL_slottime, S, numUE);
fprintf(1, '\tTransmission probability (tau_{SL})= %5.4f\n \n', prob_transmission);


