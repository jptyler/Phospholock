%%% Randomly simulate the phospholock model %%%

global m 

% Loop through all values of phosphorylation/dephosphorylation strenghs

phosphorylation_strengths = 4:-0.1:-4;

fs = 0.01;

all_parameters = zeros(10000,12);

Kds = [1 0.1, 0.01, 0.001, 0.0001];

for i = 1:10000
    
    p = zeros(12,1);
    p([1:6 end]) = 100*rand(7,1);
    p(7) = 100*rand;
    p(8) = p(7);
    p(9) = 100*rand;
    p(10) = p(9);
    p(11) = 10*rand;
    
    all_parameters(i,:) = p;
    
end

for ii = 1:length(Kds)

    all_counts = [];

    good_pars = [];

    for i = 1:length(phosphorylation_strengths)
        i
        c = 10^phosphorylation_strengths(i);
        counts = 0;

        for j = 1:10000

            % Randomize initialization
            
            initials = 100*rand(1,3);
            
            % Randomize the parameter set p

            p = all_parameters(j,:);

            p(10) = c*p(10);

            p(8) = Kds(ii)*p(7);

            [t,x] = ode23tb(@ode_phospholock, [0 100], initials, [], p);

            m = mean(x(:,1));

            options = odeset('Events', @event);

            [~,~,te,~,~] = ode23tb(@ode_phospholock, [0 1000], initials, options, p);

            if(isempty(te))
              continue;
            elseif(te(end)>900)

              if(length(te)< 5)
                  continue;
              elseif(abs((te(end)-te(end-1))-(te(end-1)-te(end-2)))<0.001)

                figure(1)

                plot(t,x(:,1))

                saveas(1, ['../Figures/phospholock_' num2str(i) '_' num2str(j) '.fig'])

                counts = counts+1;

                good_pars = [good_pars; i j];

              end

            end
        end

        all_counts = [all_counts; c counts];

    end

    save(['../Phospholock Simulation Results/Kd_eq_' num2str(Kds(ii)) '.mat'], 'good_pars', 'all_counts', 'counts', 'phosphorylation_strengths')

end

function [value, isterminal, direction] = event(~,x,~)
global m
value = x(1)-m;
isterminal = 0;
direction = 1;
end