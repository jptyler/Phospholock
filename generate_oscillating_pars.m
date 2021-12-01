% Wt = p(1); % Total WCC concentration
% a1 = p(2); % rate constant of transcription of frq
% b1 = p(3); % degradation of frq mRNA
% a2 = p(4); % translation rate of frq mRNA
% b2 = p(5); % degradation rate of FRQ
% a3 = p(6); % synthesis of FRQ complex
% b3 = p(7); % degradation of FRQ complex
% k1f = p(8); % rate of binding of unphosphorylated WCC and FCH
% k1r = p(9); % rate of unbinding of unphosphorylated WCC and FCH
% k2f = p(10); % rate of phosphorylation of unphosphorylated WCC and FCH
% k2r = p(11); % rate of dephosphorylation of phosphorylated WCC and FCH
% k3 = p(12); % rate of dissociation of phosphorylated WCC and FCH
% k4f = p(13); % rate of autophosphorylation of unphosphorylated WCC
% k4r = p(14); % rate of autodephosphorylation of phosphorylated WCC
global m

% all_pars = [];

for i = 4:10
    
    c = 1;
    
    while c <= 100

        Kd = 5;

        WT = 100*rand;
        k1f = rand;
        k1r = 10^(-Kd)*k1f;
        k2f = rand;
        k2r = .1*k2f;
        k3 = i;
        k4f = rand;
        k4r = k4f;

        p = [WT 100*rand(1,6) k1f k1r k2f k2r k3 k4f k4r];
        
        initials = rand(3,1);
        
        [t,x] = ode23tb(@ode_neuro, [0 500], initials, [], p);
        
        
        
        m = mean(x(:,1));

        options = odeset('Events', @events);

        [~,~,te,~,~] = ode23tb(@ode_neuro, [0 1000], initials, options, p);

        if(length(te)<3)

            continue;

        elseif(te(end)<900)

            continue;

        elseif(abs(te(end)-te(end-1)-te(end-1)+te(end-2))<.001)

            c = c+1

            all_pars = [all_pars; p];

            figure(1)
            plot(t,x(:,1))
        end
    end
    
    save('oscillating_params_stoichiometric_analysis.mat','all_pars')

end


function [value, isterminal, direction] = events(~,x,~)
global m
value = x(1)-m;
isterminal = 0;
direction = 1;
end

