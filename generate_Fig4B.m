%%% Plot f(R) for Kd = 1 and varying k2 ratios (k2f : k2r)

load('parameter_for_fR_plot.mat')

a1 = p(1);
b1 = p(2);
a2 = p(3);
b2 = p(4);
a3 = p(5);
b3 = p(6);
k1f = p(7);
k1r = p(8);
k2f = p(9);
k2r = p(10);
k3 = p(11);
A_T = p(12);

k2s = 0:-.01:-4;

c1 = 1:length(k2s);

c = [2*c1'/length(k2s) zeros(length(c1),1) (1-c1'/length(k2s)).^2];

c(c(:,1)>1,1) = 1;

for i = 1:length(k2s)
    
    k2r = p(10)*10^(k2s(i));

    % New Parameters

    K1 = (k2f+k2r+k3)/(k2r+k3);
    K2 = (k1r*k2r+k1r*k3+k2f*k3)/(k1f*(k2r+k3));

    % Solve for AR
    
    index = 0.01;
    
    x = 75:index:150;

    AR = (K1*A_T+K1.*x+K2-sqrt((K1*A_T+K1.*x+K2).^2-4*K1^2*A_T.*x))/(2*K1^2);

    A = A_T-K1*AR;
    
    df1 = (diff(A)./index).*(x(2:end)./A(2:end));

    figure(1)
    plot(x,A, 'Color', c(i,:))
    hold on
    plot(x(2:end),df1, 'Color', c(i,:))
    
end