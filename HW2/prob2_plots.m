% Southern Methodist University -- Math Department
% Math 6321 -- Fall 2020
% Homework 2 -- Sep 18
% Problem 2 -- forward Euler Method
% Sharon Yang -- xiny@smu.edu

% Problem Statement: Prob2
% (d) Create an overlaid plot of y1(t) and y2(t), and save this to disk.
% (e) Create a plot of y2 vs y1, and save this to disk.

% beta = 2; (d)
file2a = '2a_results.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(file2a,delimiterIn,headerlinesIn);
figure 
plot(A.data(:,1),A.data(:,2),'r');
hold on 
plot(A.data(:,1),A.data(:,3),'b');
xlabel('t');
ylabel('y');
title('Solutions over time, beta = 2');
legend('y1(t)','y2(t)');
saveas(gcf,'beta=2_(d).png');
% beta = 2; (e)
figure 
plot(A.data(:,2),A.data(:,3));
xlabel('y1');
ylabel('y2');
title('y2 vs y1, beta = 2');
saveas(gcf,'beta=2_(e).png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta = 4; (d)
file2b = '2b_results.txt';
delimiterIn = ' ';
headerlinesIn = 1;
B = importdata(file2b,delimiterIn,headerlinesIn);
figure 
plot(B.data(:,1),B.data(:,2),'r');
hold on 
plot(B.data(:,1),B.data(:,3),'b');
xlabel('t');
ylabel('y');
title('Solutions over time, beta = 4');
legend('y1(t)','y2(t)');
saveas(gcf,'beta=4_(d).png');
% beta = 4; (e)
figure 
plot(B.data(:,2),B.data(:,3));
xlabel('y1');
ylabel('y2');
title('y2 vs y1, beta = 4');
saveas(gcf,'beta=4_(e).png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta = 3.55; (d)
file2C = '2C_results.txt';
delimiterIn = ' ';
headerlinesIn = 1;
C = importdata(file2C,delimiterIn,headerlinesIn);
figure 
plot(C.data(:,1),C.data(:,2),'r');
hold on 
plot(C.data(:,1),C.data(:,3),'b');
xlabel('t');
ylabel('y');
title('Solutions over time, beta = 3.55');
legend('y1(t)','y2(t)');
saveas(gcf,'beta=3.55_(d).png');
% beta = 3.55; (e)
figure 
plot(C.data(:,2),C.data(:,3));
xlabel('y1');
ylabel('y2');
title('y2 vs y1, beta = 3.55');
saveas(gcf,'beta=3.55_(e).png');

fprintf(['When beta = 2 < 3.5(critical beta), \n'...
    'the trajectories oscillate without damping \n'...
    'and the phase space shows a stable limit cycle. \n'...
    'When beta = 4 > 3.5, '...
    'the trajectories decay in amplitude \n'...
    'and spiral in phase space into a stable fixed point. \n'...
    'When beta = 3.55 close to 3.5, \n'...
    'the amplitude of oscillation decreases as time passes \n'...
    'and approaches to a limit. '...
    'The spiral in phase space \n eventually turns into a stable cycle. '])
