function [cp_peak, cp_int, cf_peak, cf_int, ind_int] = measure_summation(filename)
%MEASURE_CABLE_SIMULATIONS allows for measurement of individual responses and
%temporal summation of responses. Measures peak amplitude and integral of
%the responses. Note that it really only works with simulated NEURON data
%generated with the cable simulations for direction sensitivity

% INPUT VARIABLES:  filename is a timeseries of voltage simulated
%                   simulated with the direction sensitivity protocal 
%                   in NEURON. Thus, it is a 12-second trace with
%                   individual responses (in first 6 seconds) and summation
%                   summation respsonses in centripetal and then
%                   centrifugal directions
%                   All files have a dt of 0.025 ms
%
% OUTPUT VARIABLES:  cp_peak = peak of response amplitude (mV; abs value)
%                               for summation in centripetal direction
%                    cp_int = integral of response (mV*s; abs value) for
%                               for summation in centripetal direction
%                    cf_peak = peak of response (mV*s; abs value) for
%                               for summation in centrifugal direction
%                    cf_int = integral of response (mV*s; abs value) for
%                               summation in centrifugal direction
%                    
%
% NOTE ON OUTPUT STRUCTURE: each output variable is a 6-element vector that
%                           corresponds to one model cable with an axial 
%                           resistance (Ra; ohm*cm) of:
%                           1, 10 50 100, 150, 200 ohm*cm, respectively
%


data = load([filename]);
dt = 0.1; % simulation time step. This is important for converting time scale into seconds

% t_vect = [0:dt:72000.125]'; 
% figure(1)
% plot(t_vect, data)

% Separate and Plot Simulated Data or 6 cables varying in axial resistance
% (Ra)
t_vect = [0:dt:12000]/1000; % in seconds
d(:,1) = data(1:length(t_vect));
d(:,2) = data(length(t_vect)+1:2*length(t_vect));
d(:,3) = data(2*length(t_vect)+1:3*length(t_vect)); 
d(:,4) = data(3*length(t_vect)+1:4*length(t_vect)); 
d(:,5) = data(4*length(t_vect)+1:5*length(t_vect)); 
%d(:,6) = data(5*length(t_vect)+1:6*length(t_vect));
%legend('Ra = 1 ohm*cm','Ra = 10 ohm*cm','Ra = 50 ohm*cm','Ra = 100 ohm*cm','Ra = 150 ohm*cm','Ra = 200 ohm*cm')

% Calculate centripetal and centrifugal integrals and peak amplitudes

d_offset = [];

for i = 1:5
    cp_peak(i) = -50-min(d(65000:90000,i));
    d_offset(:,i) = d(:,i)-(-50);
    ind_int(i) = -(trapz(d_offset(1:60000,i))*dt)/1000; % mV*s 
    cp_int(i) = -(trapz(d_offset(65000:90000,i))*dt)/1000; % mV*s 
    cf_peak(i) = -50-min(d(95000:end,i));
    cf_int(i) = -(trapz(d_offset(95000:end,i))*dt)/1000; % mV*s
end


