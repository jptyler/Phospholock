
% Load the Kd values

Kds = [1 0.1, 0.01, 0.001, 0.0001];

% iterate through all of the Kd values

for i = 1:length(Kds)
    
    load(['../Phospholock Simulation Results/Kd_eq_' num2str(Kds(i)) '.mat']) % Load each results file
    
    figure(1)
    hold on
    plot(phosphorylation_strengths, all_counts(:,2)/100)
end

figure(1)
hold off