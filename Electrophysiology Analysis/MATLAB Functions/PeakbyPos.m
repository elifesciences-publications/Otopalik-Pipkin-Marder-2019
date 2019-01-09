function [E_rev, Max_peak] = PeakbyPos(baseline_Vm_vect, peaks_vect, Vm)
%PeaksbyPos fits the peak vs. baseline Vm function with a linear regression 
%and calculates the approximate the theoretical response amplitude at a given Vm

% note that each row is a position, each column is at a a different baseline Vm 
x = baseline_Vm_vect;
y = peaks_vect;

% Linear Regression
[R,p] = corrcoef(x,y);
C = cov(x,y); slope = R(1,2)*sqrt(C(2,2)/C(1,1)); b = mean(y) - slope*mean(x);
E_rev = -b/slope;
X = [-100 -30]; % chosen to match the x-axis of interest
Y = slope * X + b;
%plot(X,Y,'k-', 'LineWidth', 1);
Y_fit = slope * x + b;
Max_peak = slope * (Vm) + b;


end

