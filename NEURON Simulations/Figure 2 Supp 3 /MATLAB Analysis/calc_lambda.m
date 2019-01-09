function [lambda_rr, resp_peak_rr, d] = calc_lambda(Filename, dir_offset)


mydir = dir(Filename);
dist_vect = [0; 150; 300; 450; 600; 750; 850];
lambda = [];
dt = 0.45; % in msec
t_vect = [0:dt:2000]/1000; % in seconds


resp_peak = zeros(7,length(mydir)-dir_offset);
lambda = zeros(1, length(mydir)-dir_offset);
for i = 1:length(mydir)-dir_offset
    clear data; clear d;
    filename = mydir(i+dir_offset).name;
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
        resp_peak(j,i) = -50-min(d(2000:3000,j)); % mV
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

% rearrange lambda and resp_peak matrices so rows correspond w/increasing
% Ra values from 50 to 300 ohm*cm:

l = lambda;
lambda_rr = [l(5); l(1:4)']; 
rp = resp_peak; 
resp_peak_rr = zeros(7,5);
resp_peak_rr(:,1) = rp(:,5);
resp_peak_rr(:,2:5) = rp(:,1:4);

end

