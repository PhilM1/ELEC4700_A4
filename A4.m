clear;
close all;
clc;
%================================
%-------- Question 1 & 2 --------
%================================
load('Carr.mat');%current array from modified A3 simulation sweep
load('Varr.mat');%voltage array from modified A3 simulation sweep

%lets get ourselves a fit to the data
P = polyfit(voltageArr, avgCurrArr, 1); %this gives us mx+b where m = P(1), b = P(2)
fittedLine = P(1)*voltageArr+P(2); %fitted line that we can plot
R3 = 1/P(1);

%Plot the current over voltage
figure(11);
plot(voltageArr,avgCurrArr);
hold on;
plot(voltageArr,fittedLine,'r');
title(['Current vs Voltage (R3 = ', num2str(R3), ' ohm)']);
xlabel('Voltage (V)');
ylabel('Current (A)');