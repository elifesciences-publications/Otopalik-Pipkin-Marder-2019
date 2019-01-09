function [baseline_Vm_mat, peaks_mat] = ErevByPos_v2(Erev_trace_mat, num_pulses, t_start, t_end)
%MEAM_RESPONSES allows for measurement of individual responses, plots as
%a function of baseline membrane potential, and calculates Erev for
%num_pulses # of positions (that are stimulated in a train at 0.5 Hz)

% INPUT VARIABLES:  filename: .abf file as string (ex. 'AGO_009__0009.abf')
%                   laser_ind is a vector of laser TTL indices from index a
%                   to b as [a:b]
%                   t_end denotes the length of the traces
% OUTPUT VARIABLES: traces: i x 26001 array of responses (i = # pulses)
%                   resp_peaks: i x 1 vector of peak amplitudes (mV)
%                   resp_int: i x 1 vector of response integrals (mV*ms)
%                   Vm_resp: lists the Vm at time of stimulation (i.e. for
%                   Erev calculation)

% calculate peak amplitude of responses:

[m,n] = size(Erev_trace_mat);

%plot ex traces each position


% detect response amplitudes for all stimuli and positions
for i = 1:m
    baseline_Vm_vect(i) = mean(Erev_trace_mat(i,4800:5000));
    peaks_vect(i) =  mean(Erev_trace_mat(i,t_start:t_end))- baseline_Vm_vect(i); % 5200:5800 is typically good time window  
end

% create matrices of baseline Vm and peak amplitudes
% where each row is for a different position.

for i = 1:num_pulses
    baseline_Vm_mat(i,:) = baseline_Vm_vect(i:num_pulses:m);
    peaks_mat(i,:) = peaks_vect(i:num_pulses:m);  
end    


end

   
    



