%MNPA for ELEC4700
%Philippe Masson
%clear;
close all;

R1 = 1;
Cap = 0.25;
R2 = 2;
L = 0.2;
R3 = 25.356630807098327; %from A3 simulation
a = 100;
R4 = 0.1;
R0 = 1000;
Vin = -10;
freq = 100;

%    I1, I2, I3,   V3,   V2,    V1,  V
G =[ -1,  1,  1,    0,    0,     0,  0;
      0,  0,  0,   -1,    1,     0,  0; 
     -1,  0,  0,    0, 1/R1, -1/R1,  0;
      0,  0, -1, 1/R3,    0,     0,  0;
      0, -1,  0,    0, 1/R2,     0,  0;     
      0,  0,  a,    0,    0,     0, -1; 
      0,  0,  0,    0,    0,     1,  0];
  
C = [ 0,  0,  0,    0,    0,     0,  0;
      0,  0,  L,    0,    0,     0,  0; 
      0,  0,  0,    0,  Cap,  -Cap,  0;
      0,  0,  0,    0,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0;    
      0,  0,  0,    0,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0];
  
F = [0;
     0;
     0;
     0;
     0;
     0;
     Vin];

%DC Sweep
Vin = -10;
V3Array = [];
V0Array = [];
VinArray = [];
Vcount = 1;
while (Vin <= 10)
    F = [0;
    0;
    0;
    0;
    0;
    0;
    Vin];
    VDC = G\F;
    V3Array(Vcount) = VDC(4);
    V0Array(Vcount) = VDC(7)*R0/(R4+R0);
    VinArray(Vcount) = Vin;
    Vcount = Vcount + 1;
    Vin = Vin + 1;
end

%plot
figure();
plot(VinArray, V0Array);
title('V0 (DC)');
xlabel('Vin (V)');
ylabel('V0 (V)');
figure();
plot(VinArray, V3Array);
title('V3 (DC)');
xlabel('Vin (V)');
ylabel('V3 (V)');


%AC Sweep;
VAC = (G+1i*2*pi*freq*C)\F;
freq = 1;
Vin = 10;
V0ACArr = [];
freqArr = [];
gainArr = [];
Vcount = 1;
 while (freq <= 1e2)
    VAC = (G+1i*2*pi*freq*C)\F;
    V0ACArr(Vcount) = VAC(7)*R0/(R4+R0);
    freqArr(Vcount) = freq;
    gainArr(Vcount) = 20*log10(V0ACArr(Vcount)/Vin);
    freq = freq + 1;
    Vcount = Vcount + 1;
 end
 
 %plot  
 figure();
 plot(freqArr, V0ACArr);
 title('V0 (AC)');
 xlabel('Frequency (Hz)');
 ylabel('V0 (V)');
 
 figure();
 plot(freqArr, gainArr);
 title('Gain (AC, 20dB per Decade)');
 xlabel('Frequency (Hz)');
 ylabel('Gain (dB)');
 
 %AC Sweep2;
 CArr = normrnd(Cap,0.05,[100000,1]);
 CArr = sort(CArr);
 
 Vin = 10;
 V0AC2Arr = [];
 gain2Arr = [];
 for c = 1:length(CArr)
    VAC = (G+1i*pi*CArr(c))\F;     
    V0AC2Arr(c) = VAC(7)*R0/(R4+R0);
    gain2Arr(c) = 20*log10(V0AC2Arr(c)/Vin);
 end
 
 %plot
 figure();
 plot(CArr, gain2Arr);
 title('Gain of C Perturbations (AC, 20dB per Decade)');
 xlabel('C (F)');
 ylabel('Gain (dB)');
 
 figure();
 histogram(real(gain2Arr))
 xlim([15.7519, 15.755]);
 title('Gain Histogram of C Perturbations');
 xlabel('Gain (dB)');
 ylabel('Count');
 
% ASSIGNMENT TRANSIENT QUESTION
%================================
%---------- Question 4 ----------
%================================
%Transient analysis
simLength = 1; %simulation length in seconds
numSteps = 1000; %number of total steps in the simulation
timestep = simLength / numSteps;
curTime = 0;
V0TransArr = [];
VinTransArr = [];
timeArr = [];
gainArr = [];
Vcount = 1;
VTransPrev = 0;
Vin2Freq = 1/0.03;
PulseOffsetTime = 0.06;
PulseStdDev = 0.03;
 while (curTime <= simLength)
    if(curTime < 0.03) %First component of input, step function
        Vin1 = 0;
    else
        Vin1 = 1;
    end        
    Vin2 = sin(2*pi*curTime*Vin2Freq); %Second component of input, sine wave   
    if(curTime >= PulseOffsetTime) %Third component of input, Gaussian Pulse
        Vin3 = exp((-(curTime - PulseOffsetTime)^2)/(2*(PulseStdDev^2)));
    else
        Vin3 = 0;
    end    
    VinNew = Vin1+Vin2+Vin3; %Combined inputs
     
    F = [0;
    0;
    0;
    0;
    0;
    0;
    VinNew];

    Fnew = (F - C*VTransPrev);
    VTrans = (G + C)\Fnew;
    VTransPrev = VTrans;
    V0TransArr(Vcount) = VTrans(7)*R0/(R4+R0);
    VinTransArr(Vcount) = VinNew;
    timeArr(Vcount) = curTime;
    gainArr(Vcount) = 20*log10(V0TransArr(Vcount)/VinNew);
    Vcount = Vcount + 1;
    curTime = curTime + timestep;
 end

%Transient Plots
figure();
plot(timeArr, V0TransArr);
hold on;
plot(timeArr, VinTransArr);
title('Vin and Vout over Time for Transient Simulation');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Vout', 'Vin');

%Frequency domain (Fourier)
f = (1/timestep)*(-(numSteps/2):(numSteps/2-1))./numSteps;
VinFourier = abs(fftshift(fft(VinTransArr)));
VoutFourier = abs(fftshift(fft(V0TransArr)));
figure();
plot(f, VoutFourier);
hold on;
plot(f, VinFourier);
title('Input and Output Frequency Response of Transient Simulation');
xlabel('Frequency (Hz)');
ylabel('FFT(Vin)');
legend('Output', 'Input');
 


% ASSIGNMENT NOISE QUESTION
%================================
%---------- Question 5 ----------
%================================
%Transient analysis with Noise

Cn = 0.00001;
NoiseMag = 0.001;

%    I1, I2, I3,   V3,   V2,    V1,  V
C = [ 0,  0,  0,    0,    0,     0,  0;
      0,  0,  L,    0,    0,     0,  0; 
      0,  0,  0,    0,  Cap,  -Cap,  0;
      0,  0,  0,   Cn,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0;    
      0,  0,  0,    0,    0,     0,  0;
      0,  0,  0,    0,    0,     0,  0];
 
simLength = 1; %simulation length in seconds
numSteps = 1000; %number of total steps in the simulation
timestep = simLength / numSteps;
curTime = 0;
V0TransArr = [];
VinTransArr = [];
timeArr = [];
gainArr = [];
Vcount = 1;
VTransPrev = 0;
Vin2Freq = 1/0.03;
PulseOffsetTime = 0.06;
PulseStdDev = 0.03;

while (curTime <= simLength)
    if(curTime >= PulseOffsetTime) %Third component of input, Gaussian Pulse
        Vin3 = exp((-(curTime - PulseOffsetTime)^2)/(2*(PulseStdDev^2)));
    else
        Vin3 = 0;
    end    
    VinNew = Vin3; %Combined inputs
    In = 0.001*normrnd(0,1);

    F = [0;
         0;
         0;
         In;
         0;
         0;
         VinNew];

    Fnew = (F - C*VTransPrev);
    VTrans = (G + C)\Fnew;
    VTransPrev = VTrans;
    V0TransArr(Vcount) = VTrans(7)*R0/(R4+R0);
    VinTransArr(Vcount) = VinNew;
    timeArr(Vcount) = curTime;
    gainArr(Vcount) = 20*log10(V0TransArr(Vcount)/VinNew);
    Vcount = Vcount + 1;
    curTime = curTime + timestep;
end

%Transient Plots with Noise
figure();
plot(timeArr, V0TransArr);
hold on;
plot(timeArr, VinTransArr);
title(['Vin and Vout over Time for Transient Noise Simulation (Cn=', num2str(Cn),')']);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Vout', 'Vin');

%Frequency domain (Fourier)
f = (1/timestep)*(-(numSteps/2):(numSteps/2-1))./numSteps;
VinFourier = abs(fftshift(fft(VinTransArr)));
VoutFourier = abs(fftshift(fft(V0TransArr)));
figure();
plot(f, VoutFourier);
hold on;
plot(f, VinFourier);
title('Input and Output Frequency Response of Transient Noise Simulation');
xlabel('Frequency (Hz)');
ylabel('FFT(Vin)');
legend('Output', 'Input');

%Transient Plots with Noise
figure();
plot(timeArr, V0TransArr);
title(['Vout over Time for Transient Noise Simulation (Cn=', num2str(Cn),')']);
xlabel('Time (s)');
ylabel('Voltage (V)');

%Transient Plots with Noise
% figure(42);
% plot(timeArr, V0TransArr1);
% hold on;
% plot(timeArr, V0TransArr2);
% plot(timeArr, V0TransArr3);
% title('Vout over Time for Transient Noise Simulations with varying Cn');
% xlabel('Time (s)');
% ylabel('Voltage (V)');
% legend('Cn = 0.1', 'Cn = 0.01', 'Cn = 0.001');

% figure(43);
% plot(timeArr1, V0TransArr1);
% hold on;
% plot(timeArr2, V0TransArr2, 'g');
% plot(timeArr3, V0TransArr3, 'r');
% title('Vout over Time for Transient Noise Simulations with varying timestep');
% xlabel('Time (s)');
% ylabel('Voltage (V)');
% legend('timestep = 0.001', 'timestep = 0.002', 'timestep = 0.01');

