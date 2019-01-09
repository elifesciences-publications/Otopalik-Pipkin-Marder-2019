function [avg_trace, std_trace, avg_peak, std_peak, avg_int, std_int] = mean_summation_v2(traces, t_end, Vm)
%MEAM_RESPONSES allows for measurement of individual responses and
%temporal summation of respones. Measures peak amplitude and integral of
%the responses.

% INPUT VARIABLES:  filename: .abf file as string (ex. 'AGO_009__0009.abf')
%                   laser_ind is a vector of laser TTL indices from index a
%                   to b as [a:b]
%                   t_end denotes the length of the traces
%                   Vm denotes the baseline membrane potential (if TECC)
% OUTPUT VARIABLES: traces: i x 26001 array of responses (i = # pulses)
%                   resp_peaks: i x 1 vector of peak amplitudes (mV)
%                   resp_int: i x 1 vector of response integrals (mV*ms)
%                   Vm_resp: lists the Vm at time of stimulation (i.e. for
%                   Erev calculation)

% calculate peak amplitude of responses:



% calculate peak and integrals

[m,n] = size(traces);
stim_ind = n;
k = 1;
for i = 1:stim_ind
    baseline_Vm(k) = mean(traces(1:500,i));
    peaks(k) =  min(traces(5000:t_end,i))- baseline_Vm(k)
    trace_offset(:,k) = traces(1:t_end,i) - baseline_Vm(k);
    integrals(k) = trapz((1:length(trace_offset(:,k)))./10000, trace_offset(:,k));
    k = k+1;
end

avg_trace = mean(trace_offset')+Vm;
std_trace = std(trace_offset');

avg_peak = mean(peaks);
std_peak = std(peaks);
avg_int = mean(integrals');
std_int = std(integrals');

% if want to plot

figure
hold on
plot((1:length(avg_trace))./10000, avg_trace+std_trace, 'c', 'LineWidth', 1)
plot((1:length(avg_trace))./10000,avg_trace-std_trace, 'c', 'LineWidth', 1)
plot((1:length(avg_trace))./10000,avg_trace,'b', 'LineWidth', 2)
ylabel('mV'); xlabel('Time (s)'); box off; 
text(2, Vm-1, ['Avg Peak = ' num2str([avg_peak]) ' mV'])
text(2, Vm-1.5, ['Std Peak = ' num2str([std_peak]) ' mV'])
text(2, Vm-2, ['Avg Int = ' num2str([avg_int]) ' mV*sec'])
text(2, Vm-2.5, ['Std Int = ' num2str([std_int]) ' mV*sec'])
ylim([mean(avg_trace)-5 mean(avg_trace)+2])
end

   
    



