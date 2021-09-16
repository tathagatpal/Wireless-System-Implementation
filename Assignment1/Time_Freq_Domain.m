fs = 30;

freq = 100;
A = 1;
phase = 0;
t = 0:1/fs:1;

x = A*sin(2*pi*freq*t + phase);

subplot(2,1,1);
plot(t,x, 'linewidth', 1.5);
hold on;
xlabel('Time (s)'); ylabel('Amplitude (units)');
title("Sine wave, with Amplitude = " + A + " unit, frequency " ...
    + freq + "Hz, and phase = " + rad2deg(phase) + " degrees");


length_x=length(x); 
n_points = 1024;
Y  = fft(x,n_points);
f = (-n_points/2 : n_points/2-1)*fs/n_points; %Making frequency axis
Y = fftshift(Y);                            %Brings the center to zero 

subplot(2,1,2);
stem(f,abs(Y), 'linewidth', 1.5); 
title('Fourier domain representation of above sine wave');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
xlim([-2*freq 2*freq]);
