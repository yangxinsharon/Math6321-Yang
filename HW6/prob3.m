% 5.11 b) Problem Statement: Plot the delta-curves for backward Euler,
% trapezoidal, and (5.13) methods.

d = [0.9, 0.5, 0.1]; % delta
C = 'color'; 
x = [0 0]; y = [-50 50]; K = 'k'; LW = 'linewidth'; FS = 'fontsize';

% backward Euler
figure
plot(y,x,K,LW,1), hold on, plot(x,y,K)
t = linspace(0, 2*pi, 1000);
z = exp(1i*t);
r = d(1)*z-1;
s = d(1)*z;
l1 = plot(r./s,C,'b',LW,2);
r = d(2)*z-1;
s = d(2)*z;
l2 = plot(r./s,C,'r',LW,2); 
r = d(3)*z-1;
s = d(3)*z;
l3 = plot(r./s,C,'g',LW,2);
legend([l1,l2,l3],'0.9','0.5','0.1')
axis([-12 12 -12 12]), axis square, grid on
title('Backward Euler',FS,16)
saveas(gcf,'Backward.png')

% trapezoidal
figure
plot(y,x,K,LW,1), hold on, plot(x,y,K)
t = linspace(0, 2*pi, 1000);
z = exp(1i*t); 
r = d(1)*z-1;
s = (d(1)*z+1)/2;
l1 = plot(r./s,C,'b',LW,2);
r = d(2)*z-1;
s = (d(2)*z+1)/2;
l2 = plot(r./s,C,'r',LW,2);
r = d(3)*z-1;
s = (d(3)*z+1)/2;
l3 = plot(r./s,C,'g',LW,2);
legend([l1,l2,l3],'0.9','0.5','0.1')
axis([-40 0 -20 20]), axis square, grid on
title('Trapezoidal',FS,16)
saveas(gcf,'Trapezoidal.png')

% 5.13
figure
plot(y,x,K,LW,1), hold on, plot(x,y,K)
t = linspace(0, 2*pi, 1000);
z = exp(1i*t);
r = (d(1).*z).^2-d(1).*z;
s = (9*(d(1).*z).^2 + 6*d(1).*z +1)/16;
l1 = plot(r./s,C,'b',LW,2);
r = (d(2).*z).^2-d(2).*z;
s = (9*(d(2).*z).^2 + 6*d(2).*z +1)/16;
l2 = plot(r./s,C,'r',LW,2);
r = (d(3).*z).^2-d(3).*z;
s = (9*(d(3).*z).^2 + 6*d(3).*z +1)/16;
l3 = plot(r./s,C,'g',LW,2);
legend([l1,l2,l3],'0.9','0.5','0.1')
axis([-10 50 -40 40]), axis square, grid on
title('5.13',FS,16)
saveas(gcf,'5.13_Method.png')

disp('sss')