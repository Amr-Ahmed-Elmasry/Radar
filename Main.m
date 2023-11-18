% Task 1
clear; clc;
figure('Name', 'Task 1');
%1) we need to specify our parameters
samplesNumber = 18; % Samples
ts = 1/(60*10^6); % Time between each sample and the other
PRI = 1/(0.2 * 10^6); % Period between each band of pulses

%2) Calculate the amplitude (The proof of the formula is included in the
%report)
Pt = 1.5 * 10^6; %Transmit Power
T = PRI/ts;  % number of instances in one cycle including zeros
tau = ts * (samplesNumber - 1); % Sample Width
aplitude = sqrt(Pt * (tau/T));
% aplitude = sqrt(Pt);

%3) Generate one wave
y = ones(1, samplesNumber) * aplitude;
n = 0:ts:ts*samplesNumber - ts;

%4) Add the zeros
for i = samplesNumber + 1 : floor(PRI/ts)
    y(i) = 0;
    n(i) = n(i-1) + ts;
end

%5) Repeat the sequance multiple times
numOfTimes = 4;
y = repmat(y,1,numOfTimes);

for i = T+1 : numOfTimes*(T)
    n(i) = n(i-1) + ts;
end

%6) Plot
plot(n,y); ylim([0 0.1]);title('Task 1: Transmitted Signal');xlabel('Time[n] (s)');

% ======================================================================= %
% ======================================================================= %
% Task 2
disp("Printing Task 2 Data:");
%1) To calculate the unambiguous range we use pulse_train funcion
% (I got the function from the book's listing)
[dt,pav,ep,prf,ru] = pulse_train(tau, PRI, Pt);
disp("Unambigous Range(m)= ");
disp(ru * 1000);
% from R we can get delta t 
delta_t = (2*ru)/(3*10^8);

%2) Calculate the range resoution we use range_resolution funcion
% (I got the function from the book's listing)
delta_R = range_resolution(tau, 's');
disp("Range Resolution(m)= ");
disp(delta_R);

% ======================================================================= %
% ======================================================================= %
% Task 3
figure('Name', 'Task 3');
disp("Printing Task 3 Data:");
%From the radar equation we can easily calculate d

%1) Get the rest of given variables
G = 10; %Gain (db)
Pr = 0.3 /1000; % Recived power (watt)
a = 0.1; %Target cross-Section (m^2)
Ae = 0.5; % Antenna effective aperture (m^2)

%2) calculate d
d = sqrt((1/(4*pi))*sqrt((Pt*G*a*Ae)/Pr));
disp("Range(m)= ");
disp(d);

%3) calculate delta_t
tt = (2*d)/(3*10^8);
disp("Time(s)= ");
disp(tt);

%4) Get the echo signal
scale = (G*a*Ae)/(4*pi*d^2)^2; % the scale of the signal
y2 = circshift(y,floor(tt/ts)) * scale;


%5) Plot
subplot(2,1,1);
plot(n,y); ylim([0 0.1]);title('Task 3: Transmitted Signal');xlabel('Time[n] (s)');

subplot(2,1,2);
plot(n,awgn(y2,2,'measured')); ylim([0 1e-10]);title('Task 3: Recived Signal');xlabel('Time[n] (s)');

% ======================================================================= %
% ======================================================================= %
% Task 4
figure('Name', 'Task 4');
%1) Remake the single pulse module in task 1
%Generate one wave
y = ones(1, samplesNumber) * aplitude;
n = 0:ts:ts*samplesNumber - ts;

%Add the zeros
for i = samplesNumber + 1 : floor(PRI/ts)
    y(i) = 0;
    n(i) = n(i-1) + ts;
end

%2) Remake the echo in task 3
y2 = circshift(y,floor(tt/ts)) * scale;

%3) Correlate both signals
convo = xcorr(awgn(y2,2,'measured'),y);
newN = -(ts*598)/2:ts:(ts*598)/2;

%4) Plot
subplot(3,1,1);
plot(n,y);xlim([0 5e-6]); ylim([0 0.1]);title('Task 4: Transmitted Signal');xlabel('Time[n] (s)');

subplot(3,1,2);
plot(n,awgn(y2,20,'measured'));xlim([0 5e-6]); ylim([0 1e-10]);title('Task 4: Recived Signal');xlabel('Time[n] (s)');

subplot(3,1,3);
plot(newN,convo);xlim([0 5e-6]); ylim([0 1e-11]);title('Task 4: Correlation Signal');xlabel('Time[n] (s)');

% ======================================================================= %
% ======================================================================= %
% Task 5
figure('Name', 'Task 5');

%1) Repeat the sequance multiple times
numOfTimes = 4;
y = repmat(y,1,numOfTimes);

for i = T+1 : numOfTimes*(T)
    n(i) = n(i-1) + ts;
end

%2) Get the echo signal
y2 = circshift(y,floor(tt/ts)) * scale;

%3) Correlate both signals
convo = xcorr(awgn(y2,2,'measured'),y);
newN = -(ts*2398)/2:ts:(ts*2398)/2;

%4) Get the average
p1 = convo(1:T+1);
p2 = convo(T+1:2*T+1);
p3 = convo(2*T+1:3*T+1);
p4 = convo(3*T+1:4*T+1);
p5 = convo(4*T+1:5*T+1);

avg = (p1+p2+p3+p4+p5)./5;
n_Avg = 0:ts:ts*T;
%4) Plot
subplot(4,1,1);
plot(n,y);xlim([-2e-5 2e-5]); ylim([0 0.1]);title('Task 5: Transmitted Signal');xlabel('Time[n] (s)');

subplot(4,1,2);
plot(n,awgn(y2,2,'measured'));xlim([-2e-5 2e-5]); ylim([0 1e-10]);title('Task 5: Recived Signal');xlabel('Time[n] (s)');

subplot(4,1,3);
plot(newN,convo);xlim([-2e-5 2e-5]); ylim([0 2.5e-11]);title('Task 5: Correlation Signal');xlabel('Time[n] (s)');

subplot(4,1,4);
plot(n_Avg,avg);xlim([-2e-5 2e-5]); ylim([0 2.5e-11]);title('Task 5: Averaging the Correlation');xlabel('Time[n] (s)');

% ======================================================================= %
% ======================================================================= %
% Task 6
figure('Name', 'Task 6 (Bonus)');
%3) Generate one wave
y = ones(1, samplesNumber) * aplitude;
n = 0:ts:ts*samplesNumber - ts;

%4) Add the zeros
for i = samplesNumber + 1 : floor(PRI/ts)
    y(i) = 0;
    n(i) = n(i-1) + ts;
end

%5) Repeat the sequance multiple times
numOfTimes = 19;
y = repmat(y,1,numOfTimes);

for i = T+1 : numOfTimes*(T)
    n(i) = n(i-1) + ts;
end
%4) Get the echo signal
scale = (G*a*Ae)/(4*pi*d^2)^2; % the scale of the signal
y2 = circshift(y,floor(tt/ts)) * scale;
disp(numOfTimes*(T));

%5) Plot
    y2_noise = awgn(y2,2,'measured');
subplot(2,1,1);
h1 = animatedline;
subplot(2,1,2);
h2 = animatedline;
for k = 1:length(n)
    addpoints(h1,n(k),y(k));
    drawnow
    addpoints(h2,n(k),y2_noise(k));
    drawnow
end



