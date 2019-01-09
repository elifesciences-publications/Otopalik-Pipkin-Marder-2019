function [traces, laser_traces, laser_ind, Vm, laser, time] = browse_summation_v4(filename, y_min, y_max, condition, num_pulses)
%BROWSE_RESPONSES allows for visualization of entire recording file and
%Laser TTL time points. It produces plots of entire file and indvidual
%responses. 

%Load relevant recording channels from the two .abf files to be compared

raw = LoadAbf(filename);
laser = raw.data.Laser_TTL;
Vm = raw.data.Vm_1;
length(Vm)
time = ((1:length(Vm))./10000)'; %in seconds

%Plot the whole recordings:

figure
subplot(2,1,1)
plot(time, Vm)
ylim([y_min y_max])
ylabel('V_m (mV)')
title([condition])
subplot(2,1,2)
plot(time, laser)
xlabel('Time (s)')
ylabel('Laser TTL')
for i = 1:2
    subplot(2,1,i)
    box off
    %xlim([1000000 1300000])
end
hold on

% Identify 'trigger' index for each response in both conditions:

laser_on = zeros(length(laser),1);
laser_on_final = zeros(length(laser),1);
laser_ind = find(laser > 0.2); 

for i = 1:length(laser_ind)
   laser_on(laser_ind(i)) = 1;
end

for i = 1:length(laser_on)-1
    if laser_on(i) + laser_on(i+1) == 1;
        laser_on_final(i) = 1;
    else laser_on_final(i) = 0;
    end
end

laser_ind = find(laser_on_final > 0);

for i = 1:length(laser_ind)-1
    if laser_ind(i+1) - laser_ind(i) > 1500;
    laser_on_final(laser_ind(i)) = 1;
    else laser_on_final(laser_ind(i)) = 0;
    end
end

%Overlay detected responses on figure(1)
laser_ind = find(laser_on_final > 0 );
subplot(2,1,1) 
hold on
plot(laser_ind./10000, ones(length(laser_ind)).*-35, 'ro')
 
% to browse responses plotting individual responses/summation:

traces = zeros(35001, round(length(laser_ind)/num_pulses));
laser_traces = zeros(35001, round(length(laser_ind)/num_pulses));
k = 1;
for i = 1:(length(laser_ind)/num_pulses)
%     supersubplot(3,5,3,i)
%     plot(laser_ind(i)-15000:laser_ind(i)+30000, Vm(laser_ind(i)-15000:laser_ind(i)+30000), 'b', 'LineWidth', 2);
%     hold on
%     %plot(laser_ind(i)-15000:laser_ind(i)+30000, laser_ind(laser_ind(i)-15000:laser_ind(i)+30000), ones(length(laser_ind)).*-60, 'go')
%     plot(laser_ind(i)-15000:laser_ind(i)+30000, (laser(laser_ind(i)-15000:laser_ind(i)+30000).*2)-65, 'r', 'LineWidth', 2);
%     title(num2str([i]))
    traces(:,i) = (Vm(laser_ind(k)-5000:laser_ind(k)+30000));
    laser_traces(:,i) = laser(laser_ind(k) - 5000: laser_ind(k)+30000);
    k = k+num_pulses;
%     ylim([-65 -50])
%     hold on
end
% traces = traces';
% [m,n] = size(traces);
% for i = 1:n
%     resp_peaks(i) = traces(i,1) - min(traces(i,:));
% end

% resp_peaks = resp_peaks'
% %plot all in higher freq:
% figure(7)
% subplot(2,1,1)
% k = 0; hold on;
% for i = 1:5
% plot((1:length(traces))+(k*2000), traces(i,:), 'LineWidth', 2)
% k = k+1;
% end
% % actual at high freq


% 
% % Calculate Average Responses Traces And Exponential Fits
% avg_ctrl_response = mean(ctrl_traces');
% % [min_ctrl t_peak_ctrl] = (min(avg_ctrl_response));
% % t_idx_ctrl = t_peak_ctrl:t_peak_ctrl+10000;
% % [fit_ctrl,yy_ctrl] =physfit('expc',t_idx_ctrl, avg_ctrl_response(t_idx_ctrl));
% avg_DA_response = mean(DA_traces');
% % [min_DA t_peak_DA] = find(min(avg_DA_response));
% % t_idx_DA = t_peak_DA:t_peak_DA+25000;
% % [fit_DA,yy_DA] =physfit('expc',t_idx_DA, avg_DA_response(t_idx_DA));
% 


%Plot Average Traces Separately





%Plot the Average Traces with Ctrl Trace Re-Scaled for Kinetics comparison

% scale_factor = (max(avg_ctrl_response) - min(avg_ctrl_response))/(max(avg_DA_response)-min(avg_DA_response));
% offset = max(avg_ctrl_response) - max((1/scale_factor).*avg_ctrl_response)
% scaled_ctrl_response = (1/scale_factor).*avg_ctrl_response + offset
% figure(6)
% plot(scaled_ctrl_response, 'k', 'LineWidth', 3)
% hold on
% plot(avg_DA_response, 'g', 'LineWidth', 3)
% ylim([-45 -39])
% xlim([-10000 30000])
% title('Position 1 Scaled Response Traces')
% ylabel('Vm (mV)')
% xlabel('Seconds')
% legend('Control','Dopamine', 'Location', 'NorthEastOutside')
