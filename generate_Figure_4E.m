%%% Stoichiometric analysis Neurospora %%%

load('oscillating_params_stoichiometric_analysis.mat')

syms x

ratios = [];

for i=1:length(all_pars)
        
        p = all_pars(i,:);

        % Save the parameters

        k4 = p(13)/p(14);
        WT = p(1);
        a1 = p(2);
        b1 = p(3);
        a2 = p(4);
        b2 = p(5); 
        a3 = p(6); 
        b3 = p(7); 
        k1f = p(8); % rate of binding of unphosphorylated WCC and FCH
        k1r = p(9); % rate of unbinding of unphosphorylated WCC and FCH
        k2f = p(10); % rate of phosphorylation of unphosphorylated WCC and FCH
        k2r = p(11); % rate of dephosphorylation of phosphorylated WCC and FCH
        k3 = p(12); % rate of dissociation of phosphorylated WCC and FCH
        k4f = p(13); % rate of autophosphorylation of unphosphorylated WCC
        k4r = p(14);

        Kd = k1r/k1f;
        k4 = k4f/k4r;
        tK1 = 1+k2f/(k2r+k3);
        tK2 = tK1+k3*k2f/(k4r*(k2r+k3));
        tK3 = (1+k4)*(Kd + k2f*k3/(k1f*(k2r+k3)));

        % Solve for R

        AR = (WT*tK1+x*tK2+tK3-sqrt((WT*tK1+x*tK2+tK3)^2-4*tK1*tK2*WT*x))/(2*tK1*tK2);

        f_W = (WT-tK2*AR)/((1+k4));

        sols = solve(f_W == b1*b2*b3*x/(a1*a2*a3));

        s = double(sols);
        try
        ratios = [ratios; s/WT k3];
        catch
            continue;
        end
        
end

boxplot(ratios(:,1), ratios(:,2))