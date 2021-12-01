function dxdt = ode_neuro(~,x,p)

% Parameters 

WT = p(1); % Total WCC concentration
a1 = p(2); % rate constant of transcription of frq
b1 = p(3); % degradation of frq mRNA
a2 = p(4); % translation rate of frq mRNA
b2 = p(5); % degradation rate of FRQ
a3 = p(6); % synthesis of FRQ complex
b3 = p(7); % degradation of FRQ complex
k1f = p(8); % rate of binding of unphosphorylated WCC and FCH
k1r = p(9); % rate of unbinding of unphosphorylated WCC and FCH
k2f = p(10); % rate of phosphorylation of unphosphorylated WCC and FCH
k2r = p(11); % rate of dephosphorylation of phosphorylated WCC and FCH
k3 = p(12); % rate of dissociation of phosphorylated WCC and FCH
k4f = p(13); % rate of autophosphorylation of unphosphorylated WCC
k4r = p(14); % rate of autodephosphorylation of phosphorylated WCC

% Parameter functions 

Kd = k1r/k1f;
k4 = k4f/k4r;
tK1 = 1+k2f/(k2r+k3);
tK2 = tK1+k3*k2f/(k4r*(k2r+k3));
tK3 = (1+k4)*(Kd + k2f*k3/(k1f*(k2r+k3)));

% Quadratic solution

WuFu = (WT*tK1+x(3)*tK2+tK3-sqrt((WT*tK1+x(3)*tK2+tK3)^2-4*tK1*tK2*WT*x(3)))/(2*tK1*tK2);

% f function

f_W = (WT-tK2*WuFu)/((1+k4));

% ODEs 

dxdt = [1/10*(a1*f_W-b1*x(1))
        1/10*(a2*x(1)-b2*x(2))
        1/10*(a3*x(2)-b3*x(3))];
%         (a4*x(3)-b4*x(4))];
    
end