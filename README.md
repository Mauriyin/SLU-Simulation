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

The simulation models the process of UEs attempting to access the channel, handling collisions through a backoff mechanism, and calculating the probability of successful transmission. The backoff mechanism adjusts dynamically based on the history of transmissions and collisions, aiming to mitigate collisions and optimize channel usage.
### SL-U simulation
Run `slu_sim.m` for SL-U simulation results shown in Table2;
### Analysis and throughput results:
Run `analysis.m` to get the analytical model results in Table 2-4.
