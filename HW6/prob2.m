% Sharon Yang HW6
C = 'color'; c = {'b','r','g','m','y','c'};
x = [0 0]; y = [-8 8]; K = 'k'; LW = 'linewidth'; FS = 'fontsize';

plot(y,x,K,LW,1), hold on, plot(x,y,K)
t = linspace(0, 2*pi, 1000);
z = exp(1i*t); r = z.^2-1;
s = (z.^2+4.*z+1)/3; plot(r./s,C,c{3},LW,2)  
axis([-3*10^(-16) 4*10^(-16) -1.5 1.5]), axis square, grid on
title('Problem 2 Stability Region',FS,16)