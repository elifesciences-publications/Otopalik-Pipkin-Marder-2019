%%% EXAMPLE ANALYSIS SCRIPT %%%%%%

% In this example, we will extract data, measure evoked response amplitudes, reversal potentials, and summation responses for photo-activated sites
% on Branch 2 of PD 44. In my google notebook, you'll find:
%   - that the corresponding .abf files for this branch are AGO_044_0008.abf, AGO_044_0009.abf, and AGO_044_0010.abf
%   - 6 positions were photo-activated on this branch.

% Plot evoked responses at individual sites (browse), save traces in matrix traces_8 (each
% row is a membrane potential time series); Note: it is possible to produce
% a "supersubplot" of all individaul evoked traces in this script (it's
% commented out because it is a slow task). One can also view individaul
% traces by plotting each row of the traces matrix. Importantly: you must
% unzip file AGO_044_0008.abf before running this function (was too large
% to upload to github without compressing).

[traces_8] = browse_responses_v4('AGO_044_0008.abf', -80, -30, 'individual');

% Detect baseline Vm and calcuate peak amplitudes for all evoked responses
% in traces_8. Note: this matrix must be inverted to work in the following
% function:

num_pulses = 6; % # of photo-activated sites = # laser pulses in the stimulus train.
%detect local minimum or maximum in specific time window during trace:
t_0 =5200; %start of time window (in indices)
t_end = 5800' %end of time window (in indices)
[baseline_Vm, peak_amplitudes] = ErevByPos_v2(traces_8', num_pulses, 5200, 5800);

% Use ErevPlotbyPos to plot deltaV (peak amplitudes) as a function of baseline Vm for each
% position and calculate reversal potential (Erev) by fitting curve with a
% linear regression. Can play around with the Vm_range to plot and fit,
% sometimes there are apparent non-linearities at the hyperpolarized
% potentials that are a result of poor current clamp hold at these low
% potentials and/or having hit the max conductance at a less hyperpolarized
% potential:

Erev = zeros(1,num_pulses); %initialize output vector
Vm_range = [1:length(baseline_Vm) - 4]; 
for i = 1:num_pulses 
    [Erev(i)] = ErevPlotbyPos(baseline_Vm(i,Vm_range), peak_amplitudes(i,Vm_range),[-100 -30], [-5 5], figure(1), subplot(2,3,i)); hold on
end   

% Calculate response amplitude predicted by linear fit of reversal potential curve using
% PeakbyPos function (in this case, determining the max response amplitude at -40
% mV). This function will also spit out the reversal potentials for each
% site (if want to skip plotting and empirical assessment with
% ErevPlotbyPos.

for i = 1:num_pulses
    [Erev(i), MaxPeak(i)] = PeakbyPos(baseline_Vm(i,Vm_range), peak_amplitudes(i,Vm_range), -40)  
end

%% Assess Centripetal and Centrifugal Summation on this Branch:

% Load traces for evoked responses to 5 Hz pulse trains in either direction
% using browse_summation_v4 function (analogous to browse_responses_v4
% above). This function will plot the responses and notate, with red
% circles, the detected evoked responses (based on time of laser pulses)
% and save the traces in a matrix (each row is a membrane potential time series)
% in the workspace:

[traces_9] = browse_summation_v4('AGO_044_0009.abf', -80, -30, 'Centripetal -50 mV', num_pulses);
[traces_10] = browse_summation_v4('AGO_044_0010.abf', -80, -30, 'Centrifugal -50 mV', num_pulses);


% Use mean_summation_v2 function to evaluate summation responses: plot and calculate means and standard
% deviations of peak amplitudes and response integrals:

[mn_trace_9, std_trace_9, mn_peak_9, std_peak_9, mn_int_9, std_int_9] = mean_summation_v2(traces_9, 30000, -50); %centripetal -50 mV
[mn_trace_10, std_trace_10, mn_peak_10, std_peak_10, mn_int_10, std_int_10] = mean_summation_v2(traces_10, 30000, -50); %centripetal -50 mV

