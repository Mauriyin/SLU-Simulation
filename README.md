# SLU-Simulation
Simulation and analysis code for the paper "Sidelink Mode 2 in Unlicensed Bands: Basic Throughput Model &amp; Validation"
### NR-U simulation
Run `nru_sim.m` for NR-U simulation results shown in Table2. 
The MATLAB code implements a Level 2 simulation of the 5G NR Unlicensed channel access. Here's a breakdown of the implementation:
1. **Initialization**:
   - Set initial parameters like the minimum backoff window size (`W0`), maximum backoff stage (`m`), retry limit (`K`), number of UEs (User Equipments) (`numUE`), and total number of slot times for the simulation duration (`max_tt`).
   - Compute `Wmax`, the maximum backoff window size, based on `m` and `W0`.
   - Initialize counters and state variables for each UE: success (`succUE`), collision (`collUE`), current backoff window size (`W_UE`), backoff stage (`stageUE`), non-transmission count (`cnt_notrans`), and transmission count (`cnt_trans`).
   - Set the backoff counter (`bc`) to 1 for each UE and start the transmission time (`tt`) counter.

2. **Simulation Loop**:
   - Run a loop until the slot time (`tt`) reaches the maximum (`max_tt`).
   - Check if there's a transmission in the previous iteration. If so, handle the transmission:
     - If there's only one UE that transmitted, it's a successful transmission. Reset the backoff window size and stage for that UE.
     - If more than one UE transmitted, it's a collision. Increase the backoff stage for each colliding UE, update their backoff window size, and select a new backoff counter value randomly within the new window.
   - Update backoff counters (`bc`) for each UE:
     - Identify the minimum backoff counter (`slot2Trans`) and the UEs that will transmit in this slot (`transUE`).
     - Increment the non-transmission count for each UE by the amount of time until the next transmission.
     - For each UE, determine if it's transmitting or not in this slot. If transmitting, reset its backoff counter and increment its transmission count. If not, decrease its backoff counter by the amount of time until the next transmission.

3. **Results Calculation**:
   - Calculate the probability of transmission (`tau_nr`) as the ratio of the mean of transmission count to the sum of the mean of transmission count and the mean of non-transmission count for all UEs.
   - Print out the probability of transmission.

The simulation models the process of UEs attempting to access the channel, handling collisions through a backoff mechanism, and calculating the probability of successful transmission. 
### SL-U simulation
Run `slu_sim.m` for SL-U simulation results shown in Table2.  
The extension of the code for 5G SL-U introduces additional complexity to model the time-continuous effects and the gap states in SL-U access. This extension involves the incorporation of sensing slots, gap states, and more detailed tracking of backoff values and SL slot boundaries. Here's a summary of the changes and new features compared to the original NR-U code:

1. **Introduction of SL Slot and Sensing Slot Times**:
    - `SL_slottime`: Time duration of an SL-U slot, with options like 250, 500, or 1000 microseconds.
    - `sensing_slottime`: Time duration of a sensing slot, set to 9 microseconds.
    - `S`: The number of sensing slots within an SL-U slot is calculated based on `SL_slottime` and `sensing_slottime`.

2. **Gap States Accounting**:
    - Gap states represent the time intervals between transmissions and are crucial for modeling the time-continuous nature of SL-U.
    - `cnt_bs`: Counts the gap states for each UE.
    - `cnt_bs0`: Counts the gap states when the backoff counter (`bc`) is zero.
    - `cnt_bs_S`: Counts the gap state 'S' transitioning from backoff state 2.

3. **Detailed Backoff Values Tracking**:
    - `cnt_bc`: An array to count backoff values for each UE, offering a more granular view of backoff stages.
    - Tracking of backoff values and gap states at the SL slot boundaries.

4. **Transmission Time Randomization**:
    - `end_trans`: Randomly selected end point within a sensing slot, introducing variability to the transmission duration and modeling the asynchronous nature of transmissions.

5. **Enhanced Transmission and Non-Transmission Tracking**:
    - `cnt_trans_other`: Counts transmissions by other UEs but not by a particular UE.
    - Detailed tracking of transmissions (`cnt_trans`) and non-transmissions (`cnt_notrans`) at the SL slot boundaries, providing a finer resolution of channel access dynamics.

6. **Probability of Transmission Calculation Adjustments**:
    - Adjustments in the calculation of the probability of transmission to account for the time-continuous nature of SL-U access and the presence of gap states.
    - Introduction of different methods for calculating the probability of transmission, such as using counting states (`pr_trans_cntstate`) and considering the SL slot boundaries (`pr_trans2cntSLslot`).

These extensions introduce a more sophisticated model for SL-U access, capturing the intricacies of time-continuous channel access and the impact of sensing and gap states on channel access dynamics. The enhanced tracking of backoff values, gap states, and SL slot boundaries, along with the randomization of transmission times, contributes to a more accurate and detailed simulation of the 5G SL-U access mechanism.
### Analysis and throughput results:
Run `analysis.m` to get the analytical model results in Table 2-4.
