% Measure lambda for 25 test set of cables:

cd '/Users/adrianeotopalik/Desktop/NEURON-7.5/Direction Sensitivity Simulations/Newest Run/Measure Lambda/Mar 21 Cable Library/'
mydir = dir('Lambda Measurements New Shunt 10mm');
dist_vect = [0; 150; 300; 450; 600; 750; 850];
lambda = [];
dt = 0.45; % in msec
t_vect = [0:dt:2000]/1000; % in seconds

debug_vect = [];
%debug!! %% this little operation determines if there are any files that
%were cut short (missing data, will be wrong vector length)
for i = 1:length(mydir)-2
    %clear data;
    filename = mydir(i+2).name;
    data = load(filename);
    [m,n] = size(data);
    debug_vect(i,1) = m;
end
%names of files to debug:
debug_ind = find(debug_vect < 31122)
for i = 1:length(debug_ind)
    debug_names{i} = mydir(debug_ind(i)+2).name;
end

%list of filenames as processed by MATLAB:
for i = 1:720
    filename_list{i} = mydir(i+2).name;
end
filename_list = filename_list';


%redo_list = mydir(redo_ind+2);
resp_peak = zeros(7,length(mydir)-2);
lambda = zeros(1, length(mydir)-2);
for i = 1:length(mydir)-2
    clear data; clear d;
    filename = mydir(i+2).name;
    data = load(filename);
    d(:,1) = data(1:length(t_vect));
    d(:,2) = data(length(t_vect)+1:2*length(t_vect));
    d(:,3) = data(2*length(t_vect)+1:3*length(t_vect)); 
    d(:,4) = data(3*length(t_vect)+1:4*length(t_vect)); 
    d(:,5) = data(4*length(t_vect)+1:5*length(t_vect)); 
    d(:,6) = data(5*length(t_vect)+1:6*length(t_vect)); 
    d(:,7) = data(6*length(t_vect)+1:7*length(t_vect)); 
    %d_offset = [];
    for j = 1:7
        resp_peak(j,i) = -60-min(d(2000:3000,j)); % mV
        %d_offset(j,i) = d(:,j)-(-50); % offsets trace to start at 0 mV, to calculate integral
        %resp_int = -trapz(t_vect, d_offset); % mV*ms 
    end
    
    % calculate lambda for this cable
    
 
    % exponential fit of voltage response (peak) vs distance:
    clear f
    f = fit(dist_vect,resp_peak(:,i),'exp1');
    y = 0.37*max(resp_peak(:,i)); % 37% of response  @ 0 microns
    % f(x) = a*exp(b*x)
    % solve for lambda x value (distance):
    coeffvals = coeffvalues(f);
    a = coeffvals(1); b = coeffvals(2);
    %y = a*exp(b*x) % e = 2.71828
    %y/a = exp(bx)
    lambda(1,i) = [log(y/a)]/b;
end

%% Rearranging Lambda Matrices:

%note: g_pas = 0.0001 (Rm = 10000) for all cables in this set.

%rearrange: % redefine l matrix for each gpas value!
l = lambda;
lambda_s10 = []; % s5 refers to a shunt section that is 5mm long.
lambda_s10(:,1) = [l(5); l(1:4)']; %diams 0.5 to 0.5
lambda_s10(:,2) = [l(15); l(11:14)'];%diams 1 to 0.5
lambda_s10(:,3) = [l(25); l(21:24)'];%diams 5 to 0.5
lambda_s10(:,4) = [l(10); l(6:9)'];%diams 10 to 0.5
lambda_s10(:,5) = [l(20); l(16:19)'];%diams 20 to 0.5

% new matrix cols = geometry and rows = Ra values

%% Plot effective lambda across parameter space (function of Ra)
%Gpas = 1e-4 / Rm = 100, 1000, 5000, 10000, 20000 ohm*cm^2
%k = [1 6 11 16 21 26];
c = colormap(jet);'
c_vect = [1 10 20 30 40 50 60];

%% Plot all together:

% Lambda vs Ra
figure(1)
for i = 1:5
    subplot(5,1,i)
    plot(Ra_short_vect(1:5), lambda_noshunt(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    plot(Ra_short_vect(1:5), lambda_s100u(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    plot(Ra_short_vect(1:5), lambda_s300u(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1)
    plot(Ra_short_vect(1:5), lambda_s500u(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    plot(Ra_short_vect(1:5), lambda_s1(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    %plot(Ra_short_vect(1:5), lambda_s5(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    %plot(Ra_short_vect(1:5), lambda_s10(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    box off; ylim([0 3000]); xlim([50 300]); ylabel('Lambda (microns)'); xlabel('Ra (Ohms*cm)');
    
end
legend('No Shunt', 'Shunt = 100um', 'Shunt = 300um', 'Shunt = 500um', 'Shunt = 1mm') % 'Shunt = 5mm')%, 'Shunt = 10mm')

for i = 1:5
    subplot(5,1,i)
    plot([0 300], [1000 1000], 'k-', 'LineWidth', 1)
    ylim([0 5000])
end

% Lambda versus Shunt Size for varying Ra values:
figure
shunt_vect = [0 100 300 500 1000]% 5000];
test_lambda_vect = [lambda_noshunt(1,1) lambda_s100u(1,1) lambda_s300u(1,1) lambda_s500u(1,1) lambda_s1(1,1)]% lambda_s5(1,1)]
plot(shunt_vect, test_lambda_vect)
ylim([0 2000]); ylabel('Lambda (microns)'); xlabel('Shunt Size (microns)'); box off

% Response amplitudes w/Distance:
%NOTE: these are not in order, but rows are by proximal diam: 
% 0.5 10 1 20 5 microns
% Cols are in this order L to R: 100 150 200 300 50
figure(2)
for i = 1:25
    subplot(5,5,i) %in plot: rows = geometry, columns = 
    plot(dist_vect, resp_peak_noshunt(:,i),'.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    plot(dist_vect, resp_peak_100um(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    plot(dist_vect, resp_peak_300um(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    plot(dist_vect, resp_peak_500um(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    plot(dist_vect, resp_peak_1mm(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    %plot(dist_vect, resp_peak_5mm(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 3);
    %plot(dist_vect, resp_peak_10mm(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1);
    ylim([0 15]); ylabel('Response (mV)'); xlabel('Distance (um)')
    box off; hold on;
end
legend('No Shunt', 'Shunt = 100um', 'Shunt = 300um', 'Shunt = 500um', 'Shunt = 1mm')%, 'Shunt = 5mm', 'Shunt = 10mm')

% All Geometries for given shunt on same axes:
figure(3)
subplot(2,1,1)
plot(Ra_short_vect(1:5), lambda_s1(:,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
plot([0 300], [1000 1000], 'k-', 'LineWidth', 1)
box off; ylim([0 3000]); xlim([50 200]); ylabel('Lambda (microns)'); xlabel('Ra (Ohms*cm)');

subplot(2,1,2)
plot(Ra_short_vect(1:5), lambda_s500u(:,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
plot([0 300], [1000 1000], 'k-', 'LineWidth', 1)
box off; ylim([0 3000]); xlim([50 200]); ylabel('Lambda (microns)'); xlabel('Ra (Ohms*cm)');
legend('0.5', '1', '5', '10', '20')




%% 100 um Shunt Section

figure
for i = 1:5
    subplot(5,1,i)
    %plot([0 300], [1000 1000], 'k-', 'LineWidth', 1); hold on
    plot(Ra_short_vect(1:5), lambda_s100u(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    ylim([0 5000]); xlim([50 300]); ylabel('Lambda (microns)'); xlabel('Ra (Ohms*cm)');
    %t = num2str([diams20(i,1) diams20(i,2)])
    %title(t); 
    box off; 
    hold on
end
%legend('lambda = 1000', ['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])]) 
subplot(5,1,1)
title('100um Shunt')

% Also plot response peaks as a function of distance:
%NOTE: these are not in order, but rows are by proximal diam: 
% 0.5 10 1 20 5 microns
figure(1)
for i = 1:30
    subplot(5,5,i) %in plot: rows = geometry, columns = 
    plot(dist_vect, resp_peak_100um(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1)
    ylim([0 20]); ylabel('Response (mV)'); xlabel('Distance (um)')
    box off; hold on;
end


%% 1mm Shunt Section

figure
for i = 1:5
    subplot(5,1,i)
    %plot([0 300], [1000 1000], 'k-', 'LineWidth', 1); hold on
    plot(Ra_short_vect(1:5), lambda_s1(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    ylim([0 5000]); xlim([50 300]); ylabel('Lambda (microns)'); xlabel('Ra (Ohms*cm)');
    %t = num2str([diams20(i,1) diams20(i,2)])
    %title(t); 
    box off; 
    hold on
end
%legend('lambda = 1000', ['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])]) 
subplot(5,1,1)
title('1mm Shunt')

% Also plot response peaks as a function of distance:
%NOTE: these are not in order, but rows are by proximal diam: 
% 0.5 10 1 20 5 microns
figure(1)
for i = 1:30
    subplot(5,5,i) %in plot: rows = geometry, columns = 
    plot(dist_vect, resp_peak_1mm(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1)
    ylim([0 20]); ylabel('Response (mV)'); xlabel('Distance (um)')
    box off; hold on;
end


%% 5mm Shunt Section

figure
for i = 1:5
    subplot(5,1,i)
    %plot([0 300], [1000 1000], 'k-', 'LineWidth', 1); hold on
    plot(Ra_short_vect(1:5), lambda_s5(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    ylim([0 5000]); xlim([50 300]); ylabel('Lambda (microns)'); xlabel('Ra (Ohms*cm)');
    %t = num2str([diams20(i,1) diams20(i,2)])
    %title(t); 
    box off; hold on;
end
%legend('lambda = 1000', ['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])]) 
subplot(5,1,1)
title('5mm Shunt')

% Also plot response peaks as a function of distance:
figure(1)
for i = 1:30
    subplot(5,5,i) %in plot: rows = geometry, columns = 
    plot(dist_vect, resp_peak_5mm(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 3)
    ylim([0 20]); ylabel('Response (mV)'); xlabel('Distance (um)')
    box off;
    hold on
end


%% 10mm Shunt Section

figure
for i = 1:5
    subplot(5,1,i)
    %plot([0 300], [1000 1000], 'k-', 'LineWidth', 1); hold on
    plot(Ra_short_vect(1:5), lambda_s10(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    ylim([0 5000]); xlim([50 300]); ylabel('Lambda (microns)'); xlabel('Ra (Ohms*cm)');
    %t = num2str([diams20(i,1) diams20(i,2)])
    %title(t); 
    box off; hold on
end
%legend('lambda = 1000', ['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])]) 
subplot(5,1,1)
title('10mm Shunt')

% Also plot response peaks as a function of distance:
%NOTE: these are not in order, but rows are by proximal diam: 
% 0.5 10 1 20 5 microns
figure(1)
for i = 1:30
    subplot(5,5,i) %in plot: rows = geometry, columns = 
    plot(dist_vect, resp_peak_10mm(:,i), '.-', 'MarkerSize', 10, 'LineWidth', 1)
    ylim([0 20]); ylabel('Response (mV)'); xlabel('Distance (um)')
    box off; hold on;
end
