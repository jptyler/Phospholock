%%% Plot the repression function and its sensitivity %%%

close all

% Parameters

AT = 8;
Kds = -2:-.01:-4;

c1 = 1:length(Kds);

c = [c1'/length(Kds) zeros(length(c1),1) (1-c1'/length(Kds)).^2];

c(c(:,1)>1,1) = 1;

index = 0.01;

for i = 1:length(Kds)
    
    Kd = 10^(Kds(i));

    x = 7:index:11;

    A = (AT-x-Kd+sqrt((AT-x-Kd).^2+4*Kd*AT))./2;
    
    df1 = (diff(A)./index).*(x(2:end)./A(2:end));
    
    figure(1)
    plot(x,A, 'Color', c(i,:))
    hold on
    
    figure(2)
    plot(x(2:end),df1, 'Color', c(i,:))
    hold on
    
end
    
    