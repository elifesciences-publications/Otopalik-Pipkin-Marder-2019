function [E_rev] = ErevPlotbyPos(baseline_Vm_vect, peaks_vect, x_lim, y_lim, fignum, subfignum)
%EREVPLOTBYPOS scatter plots response amplitudes as a function of baseline Vm, 
%fits these data with a linear regression, and calculates the approximate
%reversal potential at that position. Any data points to be excluded
%should be defined in the input variables.x_lim and y_lim are 2-element
%vectors defining the axes.




fignum
subfignum
x = baseline_Vm_vect;
y = peaks_vect;
plot([-100 -20], [0 0], 'k-'); hold on
scatter(x, y, 80, 'filled','MarkerEdgeColor', 'k', 'LineWidth',1); hold on
xlabel('V_m (mV)')
ylabel('\Delta V (mV)')
box off    

% Linear Regression
[R,p] = corrcoef(x,y);
C = cov(x,y); slope = R(1,2)*sqrt(C(2,2)/C(1,1)); b = mean(y) - slope*mean(x);
E_rev = -b/slope;
X = [-100 -30]; % chosen to match the x-axis of interest
Y = slope * X + b;
plot(X,Y,'k-', 'LineWidth', 1);
Y_fit = slope * x + b;
error = y-Y_fit;
squarederror = error.*error;% OR, error.^2 is el-by-el raising to 2nd power
meanerr = mean(squarederror);
% text(-65, 0.8, 'MSE = '); text(-55, 0.8, num2str([meanerr]))
% text(-65, 0.6, 'R = '); text(-55, 0.6, num2str([R(2,1)]))
% text(-65, 0.4, 'p = '); text(-55, 0.4, num2str([p(2,1)]))
% text(-65, 0.2, 'n = '); text(-55, 0.2, num2str([length(x)]))
text(-65, 0.5, 'Erev= '); text(-55, 0.5, num2str([E_rev]))
xlim(x_lim); ylim(y_lim)
hold on;
   


%save(file_path)


end


