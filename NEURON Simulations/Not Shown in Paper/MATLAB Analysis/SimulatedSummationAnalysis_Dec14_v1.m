%cd '/Users/adrianegotopalik/Documents/Data/Marder Lab/Newest Run'
%mydir = dir('Direction Sensitivity output_v2')
%mydir = dir('Direction Sensitivity output_v3 Second Run')

cd '/Users/adrianeotopalik/Desktop/Dec 13 Cable Library'
mydir = dir('Output')

for i = 1:length(mydir)-3
    filename = mydir(i+3).name;
    [cp_peak(i,:), cp_int(i,:), cf_peak(i,:), cf_int(i,:), ind_int(i,:)] = measure_summation(filename);
end

%debug!!
for i = 1:length(mydir)-3
    clear data;
    filename = mydir(i+3).name;
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
for i = 1:28
    filename_list{i} = mydir(i+3).name;
end
filename_list = filename_list';



% Directional Biases & Linearity
dir_bias_peaks = cp_peak./cf_peak;
dir_bias_int = cp_int./cf_int;
dir_bias_int_sub = cp_int-cf_int;
linearity_int = cp_int-ind_int;


%set up basic plotting stuff
c = colormap(parula);
c_vect = [1 10 20 30 40 50 60];
Ra_vect = [10 50 100 150 200 300];

%% PLOT DIRECTIONAL BIAS
% as a function of Ra:
% 1 plot for each gpas value, with all morphologies in diff't colors.

figure
subplot(4,1,1)
for i = ind_5000
    plot(Ra_short_vect, dir_bias_peaks(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    %ylim([-.25 .25]); ylabel('Directional Bias'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300])
    box off
end

subplot(4,1,2)
for i = ind_10000
    plot(Ra_short_vect, dir_bias_peaks(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    %ylim([-.25 .25]); ylabel('Directional Bias'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300])
    box off
end

subplot(4,1,3)
for i = ind_15000
    plot(Ra_short_vect, dir_bias_peaks(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    %ylim([-.25 .25]); ylabel('Directional Bias'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300])
    box off
end

subplot(4,1,4)
for i = ind_20000
    plot(Ra_short_vect, dir_bias_peaks(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    ylim([-.25 .25]); ylabel('Directional Bias'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300]); 
    box off
end

%%  PLOT LINEARITY
% As a function of Ra:


figure
subplot(4,1,1)
for i = ind_5000
    plot(Ra_short_vect, linearity_int(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    ylim([-15 1]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300])
    box off
end

subplot(4,1,2)
for i = ind_10000
    plot(Ra_short_vect, linearity_int(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    ylim([-15 1]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300])
    box off
end

subplot(4,1,3)
for i = ind_15000
    plot(Ra_short_vect, linearity_int(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    ylim([-15 1]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300])
    box off
end

subplot(4,1,4)
for i = ind_20000
    plot(Ra_short_vect, linearity_int(i,:), '.-', 'MarkerSize', 10, 'LineWidth', 1); hold on
    ylim([-15 1]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); 
    xlim([50 300]); 
    box off
end




figure
for i = 1:20
    subplot(5,4,i)
    plot([0 300], [0 0], 'k-', 'LineWidth', 1); hold on
    plot(Ra_vect, linearity_int([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],:)', '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    ylim([-15 1]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); xlim([50 300])
    t = num2str([diameters(i,1) diameters(i,2)])
    title(t)
    box off
end
%legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 
legend(['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])],['Rm = ' num2str([Rm_vect(6)])]) 

figure
for i = 1:20
    subplot(5,4,i)
    plot(Ra_vect, linearity_int([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],:)', '.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    ylim([-15 .5]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); plot([0 300], [0 0], 'k-', 'LineWidth', 1)
    t = num2str([diameters(i,1) diameters(i,2)])
    title(t)
end
%legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 
legend(['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])],['Rm = ' num2str([Rm_vect(6)])]) 



figure
for i = 1:20
    for j = 1:6
        subplot(5,4,i)
        plot(gpas_short_vect, linearity_int([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],j), '.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
        ylim([-15 ]); ylabel('Linearity'); plot([10^-5 10^-2], [0 0], 'k-', 'LineWidth', 1);
        xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log'); xlabel('g_p_a_s (Ohms*cm^2)'); 
        t = num2str([diameters(i,1) diameters(i,2)])
        title(t)
    end
end
legend('Ra = 10','Ra = 50','Ra = 100','Ra = 150','Ra = 200', 'Ra = 300') 
























%% Plot effective lambda across parameter space
%Gpas = 1e-4 / Rm = 100, 1000, 5000, 10000, 20000 ohm*cm^2
%k = [1 6 11 16 21 26];
k = [1 7 13 19 25 31]
figure
for i = 1:6
    Ra = i;
    subplot(6,6,k(i))
    scatter(cf_int(ind_g01,Ra), cp_int(ind_g01,Ra),30, c(c_vect(1),:), 'filled'); hold on; 
    title('Rm = 100');
    
    subplot(6,6,k(i)+1)
    scatter(cf_int(ind_g001,Ra), cp_int(ind_g001,Ra),30, c(c_vect(2),:), 'filled'); hold on; 
    title('Rm = 1000'); 
    
    subplot(6,6,k(i)+2)
    scatter(cf_int(ind_g0002,Ra), cp_int(ind_g0002,Ra),30, c(c_vect(3),:), 'filled'); hold on; 
    title('Rm = 5000'); 
    
    subplot(6,6,k(i)+3)
    scatter(cf_int(ind_g0001,Ra), cp_int(ind_g0001,Ra),30, c(c_vect(4),:), 'filled'); hold on
    title('Rm = 10000'); 
    
    subplot(6,6,k(i)+4)
    scatter(cf_int(ind_g00006,Ra), cp_int(ind_g00006,Ra),30, c(c_vect(5),:), 'filled'); hold on
    title('Rm = 15000'); 
    
    subplot(6,6,k(i)+5)
    scatter(cf_int(ind_g00005,Ra), cp_int(ind_g00005,Ra),30, c(c_vect(6),:), 'filled'); hold on
    title('Rm = 20000'); 
    
end

for i = 1:36
    subplot(6,6,i)
    plot([0 20], [0 20], 'k-', 'LineWidth', 1.5)
    xlim([0 20]); ylim([0 20]); xlabel('Integral OUT (mV)'); ylabel('Integral IN (mV)'); 
   
end

%Box plot of directional bias
figure; 
for i = 1:6
    subplot(2,3,i)
    Ra = i
    notBoxPlot(dir_bias_peaks(ind_g01,Ra), 1); hold on
    notBoxPlot(dir_bias_peaks(ind_g001,Ra), 2)
    notBoxPlot(dir_bias_peaks(ind_g0002,Ra), 3)
    notBoxPlot(dir_bias_peaks(ind_g0001,Ra), 4)
    notBoxPlot(dir_bias_peaks(ind_g00006,Ra), 5)
    notBoxPlot(dir_bias_peaks(ind_g00005,Ra), 6)
    title(['Ra = ' num2str([Ra_vect(i)])])
    ylim([1 1.3])
end
%ylim([0 6])
%% Plot directional bias across parameter space
%Gpas = 1e-4 / Rm = 100, 1000, 5000, 10000, 20000 ohm*cm^2
%k = [1 6 11 16 21 26];
k = [1 7 13 19 25 31]
figure
for i = 1:6
    Ra = i;
    subplot(6,6,k(i))
    scatter(cf_int(ind_g01,Ra), cp_int(ind_g01,Ra),30, c(c_vect(1),:), 'filled'); hold on; 
    title('Rm = 100');
    
    subplot(6,6,k(i)+1)
    scatter(cf_int(ind_g001,Ra), cp_int(ind_g001,Ra),30, c(c_vect(2),:), 'filled'); hold on; 
    title('Rm = 1000'); 
    
    subplot(6,6,k(i)+2)
    scatter(cf_int(ind_g0002,Ra), cp_int(ind_g0002,Ra),30, c(c_vect(3),:), 'filled'); hold on; 
    title('Rm = 5000'); 
    
    subplot(6,6,k(i)+3)
    scatter(cf_int(ind_g0001,Ra), cp_int(ind_g0001,Ra),30, c(c_vect(4),:), 'filled'); hold on
    title('Rm = 10000'); 
    
    subplot(6,6,k(i)+4)
    scatter(cf_int(ind_g00006,Ra), cp_int(ind_g00006,Ra),30, c(c_vect(5),:), 'filled'); hold on
    title('Rm = 15000'); 
    
    subplot(6,6,k(i)+5)
    scatter(cf_int(ind_g00005,Ra), cp_int(ind_g00005,Ra),30, c(c_vect(6),:), 'filled'); hold on
    title('Rm = 20000'); 
    
end

for i = 1:36
    subplot(6,6,i)
    plot([0 20], [0 20], 'k-', 'LineWidth', 1.5)
    xlim([0 20]); ylim([0 20]); xlabel('Integral OUT (mV)'); ylabel('Integral IN (mV)'); 
   
end

%Box plot of directional bias
figure; 
for i = 1:6
    subplot(2,3,i)
    Ra = i
    notBoxPlot(dir_bias_peaks(ind_g01,Ra), 1); hold on
    notBoxPlot(dir_bias_peaks(ind_g001,Ra), 2)
    notBoxPlot(dir_bias_peaks(ind_g0002,Ra), 3)
    notBoxPlot(dir_bias_peaks(ind_g0001,Ra), 4)
    notBoxPlot(dir_bias_peaks(ind_g00006,Ra), 5)
    notBoxPlot(dir_bias_peaks(ind_g00005,Ra), 6)
    title(['Ra = ' num2str([Ra_vect(i)])])
    ylim([1 1.3])
end
%ylim([0 6])

%% Plot linearity of summation across parameter space
% Gpas = 1e-4 / Rm = 100, 1000, 5000, 10000, 20000 ohm*cm^2
k = [1 7 13 19 25 31]
figure
for i = 1:6
    Ra = i;
    subplot(6,6,k(i))
    scatter(ind_int(ind_g01,Ra), cp_int(ind_g01,Ra),40, c(c_vect(1),:), 'filled')
    hold on; 
    title('Rm = 100')
    subplot(6,6,k(i)+1)
    scatter(ind_int(ind_g001,Ra), cp_int(ind_g001,Ra),40, c(c_vect(2),:), 'filled')
    hold on; 
    title('Rm = 1000')
    subplot(6,6,k(i)+2)
    scatter(ind_int(ind_g0002,Ra), cp_int(ind_g0002,Ra),40, c(c_vect(3),:), 'filled')
    hold on; 
    title('Rm = 5000')
    subplot(6,6,k(i)+3)
    scatter(ind_int(ind_g0001,Ra), cp_int(ind_g0001,Ra),40, c(c_vect(4),:), 'filled')
    hold on; 
    title('Rm = 10000')
    subplot(6,6,k(i)+4)
    scatter(ind_int(ind_g00006,Ra), cp_int(ind_g00006,Ra),40, c(c_vect(5),:), 'filled')
    hold on; 
    title('Rm = 15000')
    subplot(6,6,k(i)+5)
    scatter(ind_int(ind_g00005,Ra), cp_int(ind_g00005,Ra),40, c(c_vect(5),:), 'filled')
    hold on; 
    title('Rm = 20000')
end


for i = 1:36
    subplot(6,6,i)
    plot([0 30], [0 30], 'k-', 'LineWidth', 1.5)
    xlim([0 30]); ylim([0 30]); xlabel('Expected (mV)'); ylabel('Measured (mV)');
end

% Consider Directional Sensitivity and Summation for Particular
% Morphologies

%Straight, varying diameters: 0.5, 1, 5, 10


% Straight (s) neurites:
figure
for i = 1:8
    subplot(2,4,i)
    plot([-40 0], [-40 0], 'k-', 'LineWidth',2); hold on
    box off
    if i < 5
        xlabel('OUT')
        ylabel('IN')
    else xlabel('Expected'); ylabel('Measured')
    end
end
for i = 1:6
    subplot(2,4,1)
    scatter(-cf_int_new(s05_ind, i), -cp_int_new(s05_ind, i), 40, 'filled'); hold on
    subplot(2,4,2)
    scatter(-cf_int_new(s1_ind, i), -cp_int_new(s1_ind, i), 40, 'filled'); hold on
    subplot(2,4,3)
    scatter(-cf_int_new(s5_ind, i), -cp_int_new(s5_ind, i), 40, 'filled'); hold on
    subplot(2,4,4)
    scatter(-cf_int_new(s10_ind, i), -cp_int_new(s10_ind, i), 40, 'filled'); hold on
end
for i = 1:6
    subplot(2,4,5)
    scatter(-ind_int_new(s05_ind, i), -cp_int_new(s05_ind, i), 40, 'filled'); hold on
    subplot(2,4,6)
    scatter(-ind_int_new(s1_ind, i), -cp_int_new(s1_ind, i), 40, 'filled'); hold on
    subplot(2,4,7)
    scatter(-ind_int_new(s5_ind, i), -cp_int_new(s5_ind, i), 40, 'filled'); hold on
    subplot(2,4,8)
    scatter(-ind_int_new(s10_ind, i), -cp_int_new(s10_ind, i), 40, 'filled'); hold on
end

for i = 1:8
    subplot(2,4,i)
    ylim([-20 0]); xlim([-20 0])
end

% Neurites with increasing taper (tip = 0.5 constant, end increases from
% 0.5 to 10: 0.5:0.5; 0.5:1, 0.5:5, 0.5:10

figure
for i = 1:8
    subplot(2,4,i)
    plot([-40 0], [-40 0], 'k-', 'LineWidth',2); hold on
    box off
    if i < 5
        xlabel('OUT')
        ylabel('IN')
    else xlabel('Expected'); ylabel('Measured')
    end
end

for i = 1:6
    subplot(2,4,1)
    scatter(-cf_int(s05_ind, i), -cp_int(s05_ind, i), 40, 'filled'); hold on
    subplot(2,4,2)
    scatter(-cf_int(t05_1_ind, i), -cp_int(t05_1_ind, i), 40, 'filled'); hold on
    subplot(2,4,3)
    scatter(-cf_int(t05_5_ind, i), -cp_int(t05_5_ind, i), 40, 'filled'); hold on
    subplot(2,4,4)
    scatter(-cf_int(t05_10_ind, i), -cp_int(t05_10_ind, i), 40, 'filled'); hold on
end
for i = 1:6
    subplot(2,4,5)
    scatter(-ind_int(s05_ind, i), -cp_int(s05_ind, i), 40, 'filled'); hold on
    subplot(2,4,6)
    scatter(-ind_int(t05_1_ind, i), -cp_int(t05_1_ind, i), 40, 'filled'); hold on
    subplot(2,4,7)
    scatter(-ind_int(t05_5_ind, i), -cp_int(t05_5_ind, i), 40, 'filled'); hold on
    subplot(2,4,8)
    scatter(-ind_int(t05_10_ind, i), -cp_int(t05_10_ind, i), 40, 'filled'); hold on
end

for i = 1:8
    subplot(2,4,i)
    ylim([-20 0]); xlim([-20 0])
end

% consider directional bias to reduce dimensionality:

figure
for i = 1:4
    subplot(1,4,i)
    
    box off
    if i < 5
        xlabel('OUT')
        ylabel('IN')
    else xlabel('Expected'); ylabel('Measured')
    end
end

for i = 1:6
    subplot(1,4,1)
    plot(i*ones(length(dir_bias_int(s05_ind, i)),1), dir_bias_int(s05_ind, i), 'o-'); hold on
    subplot(1,4,2)
    plot(i*ones(length(dir_bias_int(t05_1_ind, i)),1), dir_bias_int(t05_1_ind, i), 'o-'); hold on
    subplot(1,4,3)
    plot(i*ones(length(dir_bias_int(t05_5_ind, i)),1), dir_bias_int(t05_5_ind, i), 'o-'); hold on
    subplot(1,4,4)
    plot(i*ones(length(dir_bias_int(t05_10_ind, i)),1), dir_bias_int(t05_10_ind, i), 'o-'); hold on
end

% Only consider 0.5:10 tapered neurites:
% function of Ra:
figure

for i = 1:6
    scatter(Ra_vect(i)*ones(length(dir_bias_int(t05_10_ind, i)),1), dir_bias_int(t05_10_ind, i), 40, 'filled'); hold on
end
legend('Ra = 10','Ra = 50','Ra = 100','Ra = 150','Ra = 200', 'Ra = 300') 


%%%%%%%% DIRECTIONAL BIAS INTEGRALS %%%%%%%%%%%%%%%%%

%%%%%%%%%% Selected Neurites %%%%%%%%%%%%%%%%%%%%%%%

% Selected Neurites:

dir_bias_int_sub = cp_int-cf_int;
linearity_int = cp_int-ind_int;

% Selected Neurites Directional Bias:
figure
%gpas_ind = [t05_10_ind]; % to call measurements for different gpas values
subplot(2,4,1)
for i = 1:6
    plot(Ra_vect, dir_bias_int_sub(s05_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(s05_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)');title('0.5:0.5')


subplot(2,4,2)
for i = 1:6
    plot(Ra_vect, dir_bias_int_sub(t05_1_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(t05_1_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)'); title('0.5:1')


subplot(2,4,3)
for i = 1:6
    plot(Ra_vect, dir_bias_int_sub(t05_5_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(t05_5_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)'); title('0.5:5')


subplot(2,4,4)
for i = 1:6
    plot(Ra_vect, dir_bias_int_sub(t05_10_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(t05_10_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)'); title('0.5:10')
legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 

%function of gpas:

subplot(2,4,5)
for i = 1:6
    plot(gpas_vect(s05_ind), dir_bias_int_sub(s05_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(s05_ind), dir_bias_int_sub(s05_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,6)
for i = 1:6
    plot(gpas_vect(t05_1_ind), dir_bias_int_sub(t05_1_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(t05_1_ind), dir_bias_int_sub(t05_1_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,7)
for i = 1:6
    plot(gpas_vect(t05_5_ind), dir_bias_int_sub(t05_5_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(t05_5_ind), dir_bias_int_sub(t05_5_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,8)
for i = 1:6
    plot(gpas_vect(t05_10_ind), dir_bias_int_sub(t05_10_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(t05_10_ind), dir_bias_int_sub(t05_10_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

legend('Ra = 10','Ra = 50','Ra = 100','Ra = 150','Ra = 200', 'Ra = 300') 

for i = 1:4
    subplot(2,4,i)
    plot([0 300], [0 0], 'k-', 'LineWidth', 1)
    ylim([-0.3 0.3]); 
    box off
end
for i = 5:8
    subplot(2,4,i)
    plot([10^-5 10^-2], [0 0], 'k-', 'LineWidth', 1)
    ylim([-0.3 0.3]); 
    box off
end

%%%% Selected Neurites: Linearity

figure
%gpas_ind = [t05_10_ind]; % to call measurements for different gpas values
subplot(2,4,1)
for i = 1:6
    plot(Ra_vect, linearity_int(s05_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(s05_ind(i),:), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('Ra (Ohms*cm)');title('0.5:0.5')


subplot(2,4,2)
for i = 1:6
    plot(Ra_vect, linearity_int(t05_1_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(t05_1_ind(i),:), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('Ra (Ohms*cm)'); title('0.5:1')


subplot(2,4,3)
for i = 1:6
    plot(Ra_vect, linearity_int(t05_5_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(t05_5_ind(i),:), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('Ra (Ohms*cm)'); title('0.5:5')


subplot(2,4,4)
for i = 1:6
    plot(Ra_vect, linearity_int(t05_10_ind(i),:),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(Ra_vect, dir_bias_int_sub(t05_10_ind(i),:), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('Ra (Ohms*cm)'); title('0.5:10')
legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 

%function of gpas:

subplot(2,4,5)
for i = 1:6
    plot(gpas_vect(s05_ind), linearity_int(s05_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(s05_ind), dir_bias_int_sub(s05_ind, i), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,6)
for i = 1:6
    plot(gpas_vect(t05_1_ind), linearity_int(t05_1_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(t05_1_ind), dir_bias_int_sub(t05_1_ind, i), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,7)
for i = 1:6
    plot(gpas_vect(t05_5_ind), linearity_int(t05_5_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(t05_5_ind), dir_bias_int_sub(t05_5_ind, i), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,8)
for i = 1:6
    plot(gpas_vect(t05_10_ind), linearity_int(t05_10_ind, i),'.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    %scatter(gpas_vect(t05_10_ind), dir_bias_int_sub(t05_10_ind, i), 40, 'filled'); hold on
end
ylabel('Linearity')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

legend('Ra = 10','Ra = 50','Ra = 100','Ra = 150','Ra = 200', 'Ra = 300') 

for i = 1:4
    subplot(2,4,i)
    plot([0 300], [0 0], 'k-', 'LineWidth', 1)
    ylim([-10 1]); 
    box off
end
for i = 5:8
    subplot(2,4,i)
    plot([10^-5 10^-2], [0 0], 'k-', 'LineWidth', 1)
    ylim([-10 1]); 
    box off
end


%%%%%%%%%%%%%% ALL NEURITE MORPHOLOGIES %%%%%%%%%%%%%%

% All Neurite Morphologies between 0.5 and 20 microns:
% as a function of Ra: 

figure
for i = 1:20
    subplot(5,4,i)
    plot([0 300], [0 0], 'k-', 'LineWidth', 1); hold on
    plot(Ra_vect, dir_bias_int_sub([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],:)', '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    ylim([-.25 .25]); ylabel('Directional Bias'); xlabel('Ra (Ohms*cm)'); 
    t = num2str([diameters(i,1) diameters(i,2)]); xlim([50 300])
    title(t)
    box off
end
%legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 
legend(['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])],['Rm = ' num2str([Rm_vect(6)])]) 

figure
for i = 1:20
    for j = 1:6
        subplot(5,4,i)
        plot(gpas_short_vect, dir_bias_int_sub([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],j), '.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
        ylim([-.3 .3]); ylabel('Directional Bias'); plot([10^-5 10^-2], [0 0], 'k-', 'LineWidth', 1);
        xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log'); xlabel('g_p_a_s (Ohms*cm^2)'); 
        t = num2str([diameters(i,1) diameters(i,2)])
        title(t)
    end
end
legend('Ra = 10','Ra = 50','Ra = 100','Ra = 150','Ra = 200', 'Ra = 300') 

% LINEARITY
figure
for i = 1:20
    subplot(5,4,i)
    plot([0 300], [0 0], 'k-', 'LineWidth', 1); hold on
    plot(Ra_vect, linearity_int([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],:)', '.-', 'MarkerSize', 10, 'LineWidth', 1); 
    ylim([-15 1]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); xlim([50 300])
    t = num2str([diameters(i,1) diameters(i,2)])
    title(t)
    box off
end
%legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 
legend(['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])],['Rm = ' num2str([Rm_vect(6)])]) 

figure
for i = 1:20
    subplot(5,4,i)
    plot(Ra_vect, linearity_int([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],:)', '.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
    ylim([-15 .5]); ylabel('Linearity'); xlabel('Ra (Ohms*cm)'); plot([0 300], [0 0], 'k-', 'LineWidth', 1)
    t = num2str([diameters(i,1) diameters(i,2)])
    title(t)
end
%legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 
legend(['Rm = ' num2str([Rm_vect(1)])],['Rm = ' num2str([Rm_vect(2)])], ['Rm = ' num2str([Rm_vect(3)])], ['Rm = ' num2str([Rm_vect(4)])],['Rm = ' num2str([Rm_vect(5)])],['Rm = ' num2str([Rm_vect(6)])]) 



figure
for i = 1:20
    for j = 1:6
        subplot(5,4,i)
        plot(gpas_short_vect, linearity_int([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1],j), '.-', 'MarkerSize', 20, 'LineWidth', 2); hold on
        ylim([-15 ]); ylabel('Linearity'); plot([10^-5 10^-2], [0 0], 'k-', 'LineWidth', 1);
        xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log'); xlabel('g_p_a_s (Ohms*cm^2)'); 
        t = num2str([diameters(i,1) diameters(i,2)])
        title(t)
    end
end
legend('Ra = 10','Ra = 50','Ra = 100','Ra = 150','Ra = 200', 'Ra = 300') 


%%%%%%%% DIRECTIONAL BIAS PEAKS %%%%%%%%%%%%%%%%%

figure
%gpas_ind = [t05_10_ind]; % to call measurements for different gpas values
subplot(2,4,1)
for i = 1:6
    scatter(Ra_vect, dir_bias_peaks(s05_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)');title('0.5:0.5')


subplot(2,4,2)
for i = 1:6
    scatter(Ra_vect, dir_bias_peaks(t05_1_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)'); title('0.5:1')

subplot(2,4,3)
for i = 1:6
    scatter(Ra_vect, dir_bias_peaks(t05_5_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)'); title('0.5:5')


subplot(2,4,4)
for i = 1:6
    scatter(Ra_vect, dir_bias_peaks(t05_10_ind(i),:), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('Ra (Ohms*cm)'); title('0.5:10')
legend('gpas = 0.00006', 'gpas = 0.00005','gpas = 0.0001' , 'gpas = 0.0002', 'gpas = 0.001','gpas = 0.01') 



%function of gpas:

subplot(2,4,5)
for i = 1:6
    scatter(gpas_vect(s05_ind), dir_bias_peaks(s05_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,6)
for i = 1:6
    scatter(gpas_vect(t05_1_ind), dir_bias_peaks(t05_1_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,7)
for i = 1:6
    scatter(gpas_vect(t05_5_ind), dir_bias_peaks(t05_5_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

subplot(2,4,8)
for i = 1:6
    scatter(gpas_vect(t05_10_ind), dir_bias_peaks(t05_10_ind, i), 40, 'filled'); hold on
end
ylabel('Directional Bias')
xlabel('g_p_a_s (Ohms*cm^2)')
xlim([10^(-4.5) 10^-2]); set(gca, 'XScale', 'log')

legend('Ra = 10','Ra = 50','Ra = 100','Ra = 150','Ra = 200', 'Ra = 300') 

for i = 1:4
    subplot(2,4,i)
    plot([0 300], [1 1], 'k-', 'LineWidth', 1.5)
    box off; ylim([0.75 1.25]); 
end
for i = 5:8
    subplot(2,4,i)
    plot([10^-5 10^-2], [1 1], 'k-', 'LineWidth', 1.5)
    ylim([0.75 1.25]); 
    box off
end

%%%%%%%%%%%%%%%%

% HeatMaps
figure
for i = 1:20 %i= [1 2 4 3 9 10 12 11 17 18 20 19 5 6 8 7 13 14 16 15]
subplot(5,4,i)
h=heatmap(Ra_vect, gpas_vect([1+i-1, 21+i-1, 41+i-1, 61+i-1, 81+i-1, 101+i-1]), cdata([1+i-1,21+i-1,41+i-1,61+i-1,81+i-1,101+i-1],:), 'ColorLimits', [-0.25 0.25], 'Colormap', cool) 
colormap(cool)
% ax = h.plot; % 'ax' will be a handle to a standard MATLAB axes.
%h.ColorLimits = [0.9 1.1];
% colorbar('Peer', ax); % Turn the colorbar on
% caxis(ax, [0.085 1.015]);
end

