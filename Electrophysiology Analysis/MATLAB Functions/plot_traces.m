function [ output_args ] = plot_traces(trace_mat, num_pulses)
%PLOT_TRACES plots all traces in the specified matrix of traces. The number
%of positions (or pulses) in the train must also be specified.

%plot all traces on supersubplot
    k=1;
    [m,n] = size(trace_mat);
    
    for i = 1:n
        supersubplot(3,5,4,i);
        plot(trace_mat(i,:), 'k-', 'LineWidth', 2)
        xlim([0.25e4 1.25e4]); box off; 
        baseline_Vm = mean(trace_mat(i,4000:4800)); %baseline right before laser pulse
        y_max = baseline_Vm+4;
        y_min = baseline_Vm-4;
        ylim([y_min y_max])
        pos_num = k;
        title(['Pos ' num2str([pos_num])])
        if k < num_pulses
            k = k+1;
        else k = 1;
        end
    end
    
end

