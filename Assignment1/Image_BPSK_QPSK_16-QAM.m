image = imread('assign1.jpg');
%Input image of size 64x64

%Converting image to bits
stream = reshape((dec2bin(image,8)-'0').',1,[]);
[m, n] = size(image);

%Shifting the MSB to left
for i = 1:8:length(stream)-8
    stream(i:i+7) = stream(i+7 : -1: i);
end

%Uncomment to plot the bitstream of the image signal
% figure(1);
% stairs(1:1:120, stream(1:120), 'linewidth', 1.5);
% hold on;
% xlabel('Time'); ylabel('Amplitude');
% title('Time domain representation of the image');
% figure(); %Buffer figure

choice = input('Enter 1 for BPSK, 2 for QPSK, 3 for 16-QAM:  ');

if choice == 1
    %BPSK MODULATION___________________________________________________________
    bpsk_stream = zeros(1,length(stream)/2);
    bpsk_mod = pskmod(stream, 2, pi);

    %Fig. 2
    scatterplot(bpsk_mod);

    %Time domain
    %carrier frequency of 1 Hz
    t = 0:0.01:1;
    fc = 1;
    time = 0:0.01:24+0.23;
    bpsk_time = [];
    for i = 1:24
        in_phase = real(bpsk_mod(i));
        quad = imag(bpsk_mod(i));
        carrier = in_phase*cos(2*pi*fc*t) - quad*sin(2*pi*fc*t);
        bpsk_time = [bpsk_time carrier];
    end

    figure(3);
    subplot(3,1,1);
    stairs(0:23, stream(1:24), 'linewidth',1.5);
    hold on;
    % plot(time,bpsk_time, 'linewidth',1.5);
    grid on; grid minor
    xlabel('Time (s)'); ylabel('Bit level');
    title('Image bits waveform');

    subplot(3,1,2);
    plot(time,bpsk_time, 'linewidth',1.5);
    hold on;
    grid on; grid minor
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Time-domain representation of signal modulated using BPSK');
    
    %Frequency domain
    fs = 100;
    length_x=length(bpsk_time); 
    n_points = 1024;
    Y  = fft(bpsk_time,n_points);
    f = (-n_points/2 : n_points/2-1)*fs/n_points;
    Y = fftshift(Y);

    subplot(3,1,3);
    stem(f,abs(Y)); 
    hold on;
    xlabel('Frequency(Hz)'); ylabel('Amplitude');
    title('Freq-domain representation of signal modulated using BPSK');
    
    %BPSK DEMODULATION_________________________________________________________
    
%For high SNR
    high_snr = 25;
    bpsk_mod = awgn(bpsk_mod, high_snr);
    bpsk_demod = pskdemod(bpsk_mod, 2, pi);
    
    %Image reconstruction
    for i = 1:8:length(bpsk_demod)-8
        bpsk_demod(i:i+7) = bpsk_demod(i+7 : -1: i);
    end

    s_bpsk = num2cell(reshape(bpsk_demod,8,[])',2);
    b_bpsk = cellfun(@(bpsk_demod) bin2dec(strrep(num2str(bpsk_demod),' ','')), s_bpsk);
    img_bpsk=reshape(b_bpsk,m,n);

    figure(4);
    subplot(2,1,1);
    imshow(img_bpsk, [0 255]);
    hold on;
    title("Reconstructed image using BPSK, SNR = " + high_snr);


    low_snr = 5;
    bpsk_mod = awgn(bpsk_mod, low_snr);
    bpsk_demod = qamdemod(bpsk_mod, 2);

    for i = 1:8:length(bpsk_demod)-8
        bpsk_demod(i:i+7) = bpsk_demod(i+7 : -1: i);
    end

    s_bpsk = num2cell(reshape(bpsk_demod,8,[])',2);
    b_bpsk = cellfun(@(bpsk_demod) bin2dec(strrep(num2str(bpsk_demod),' ','')), s_bpsk);
    img_bpsk=reshape(b_bpsk,m,n);

    figure(4);
    subplot(2,1,2);
    imshow(img_bpsk, [0 255]);
    hold on;
    title("Reconstructed image using BPSK, SNR = " + low_snr);
    
    scatterplot(bpsk_mod);
    
    

elseif choice == 2
    %QPSK MODULATION___________________________________________________________
    qpsk_stream = zeros(1, length(stream)/2);

    %combining 2 bits into 1 as 2 bits per symbol
    % are being modulated in QPSK
    for i = 1:2:length(stream)
        qpsk_stream(round(i/2)) = bin2dec(strcat(num2str(stream(i)),num2str(stream(i+1))));
    end

    %QPSK Modulation
    qpsk_mod = pskmod(qpsk_stream,4,pi/2);
    scatterplot(qpsk_mod);                      %Fig. 3

    %Time domain
    %carrier frequency of 1 Hz
    t = 0:0.001:1;
    fc = 1;
    time = 0:0.001:24+0.023;
    qpsk_time = [];
    for i = 1:24
        in_phase = real(qpsk_mod(i));
        quad = imag(qpsk_mod(i));
        carrier = in_phase*cos(2*pi*fc*t) - quad*sin(2*pi*fc*t);

        qpsk_time = [qpsk_time carrier];
    end

    figure(3);
    subplot(2,1,1);
    plot(time,qpsk_time, 'linewidth',1.5);
    hold on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Time-domain representation of signal modulated using QPSK');
    
    %Freq domain
    fs = 1000;
    length_x=length(qpsk_time); 
    n_points = 1024;
    Y  = fft(qpsk_time,n_points);
    f = (-n_points/2 : n_points/2-1)*fs/n_points;
    Y = fftshift(Y);

    subplot(2,1,2);
    stem(f,abs(Y)); 
    hold on;
    xlabel('Frequency (Hz)'); ylabel('Amplitude');
    xlim([-100*fc 100*fc]);
    title('Freq-domain representation of signal modulated using QPSK');
    
    
    %QPSK DEMODULATION________________________________________________________
    high_snr = 25;
    qpsk_mod = awgn(qpsk_mod, high_snr); %High SNR
    qpsk_demod = pskdemod(qpsk_mod, 4, pi/2);

    qpsk_reconstruct = zeros(1,length(stream));

    for i = 1:length(qpsk_mod)
        j = i*2;
        temp = de2bi(qpsk_demod(i),2,'left-msb');

        qpsk_reconstruct(j-1) = temp(1);
        qpsk_reconstruct(j) = temp(2);

    end

%reconstruction of image
    for i = 1:8:length(qpsk_reconstruct)-8
        qpsk_reconstruct(i:i+7) = qpsk_reconstruct(i+7 : -1: i);
    end

    s = num2cell(reshape(qpsk_reconstruct,8,[])',2);
    b = cellfun(@(qpsk_reconstruct) bin2dec(strrep(num2str(qpsk_reconstruct),' ','')), s);
    img_qpsk=reshape(b,m,n);

    figure(4);
    subplot(2,1,1);
    imshow(img_qpsk, [0 255]);
    hold on;
    title("Reconstructed image using QPSK, SNR = " + high_snr);


    low_snr = 5;
    qpsk_mod = awgn(qpsk_mod, low_snr); %Low SNR
    qpsk_demod = pskdemod(qpsk_mod, 4, pi/2);

    qpsk_reconstruct = zeros(1,length(stream));

    for i = 1:length(qpsk_mod)
        j = i*2;
        temp = de2bi(qpsk_demod(i),2,'left-msb');

        qpsk_reconstruct(j-1) = temp(1);
        qpsk_reconstruct(j) = temp(2);

    end


    for i = 1:8:length(qpsk_reconstruct)-8
        qpsk_reconstruct(i:i+7) = qpsk_reconstruct(i+7 : -1: i);
    end

    s = num2cell(reshape(qpsk_reconstruct,8,[])',2);
    b = cellfun(@(qpsk_reconstruct) bin2dec(strrep(num2str(qpsk_reconstruct),' ','')), s);
    img_qpsk=reshape(b,m,n);

    figure(4);
    subplot(2,1,2);
    imshow(img_qpsk, [0 255]);
    hold on;
    title("Reconstructed image using QPSK, SNR = " + low_snr);    
    
    
elseif choice == 3
    
    %16-QAM MODULATION_________________________________________________________
    qam_stream = zeros(1, length(stream)/4);

    %Combining 4 bits as 4 bits per symbol are being transmitted in 
    %16-QAM
    for i = 1:4:length(stream)
        qam_stream(ceil(i/4)) = bin2dec(strcat(num2str(stream(i)),num2str(stream(i+1)), ...
            num2str(stream(i+2)),num2str(stream(i+3))));
    end
   
    %16-QAM Modulation
    qam16=qammod(qam_stream,16);
    scatterplot(qam16);                           %Fig. 4

    
    %Time domain
    %carrier frequency of 10 Hz
    t = 0:0.001:1;
    fc = 10;
    time = 0:0.001:4+0.003;
    qam16_time = [];
    for i = 1:4
        in_phase = real(qam16(i));
        quad = imag(qam16(i));
        carrier = in_phase*cos(2*pi*fc*t) - quad*sin(2*pi*fc*t);

        qam16_time = [qam16_time carrier/sqrt(2)];
        qam16_time(1:8) = qam16_time(1:8)/1.85;
    end
    qam16_time(2997:3002) = qam16_time(2997:3002)/1.85;
   

    figure(3);
    qam16_time(1002:2002) = qam16_time(1002:2002)/1.85;
    subplot(2,1,1);
    qam16_time(2997:4004) = qam16_time(2997:4004)/1.85;
    plot(time,qam16_time, 'linewidth',1.5);
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Time-domain representation of 16-QAM');
    
    fs = 1000;
    length_x=length(qam16_time); 
    n_points = 1024;
    Y  = fft(qam16_time,n_points);
    f = (-n_points/2 : n_points/2-1)*fs/n_points;
    Y = fftshift(Y);

    subplot(2,1,2);
    stem(f,abs(Y)); 
    hold on;
    xlabel('Frequency (Hz)'); ylabel('Amplitude');
    title('Freq-domain representation of signal modulated using QPSK');
    
    
    %16-QAM DEMODULATION_______________________________________________________

    high_snr = 25;
    qam16 = awgn(qam16, high_snr); %High SNR
    qam16_demod = qamdemod(qam16, 16);

    qam16_reconstruct = zeros(1,length(stream));

    for i = 1:length(qam16_demod)
        j = i*4;
        temp = de2bi(qam16_demod(i),4,'left-msb');

        qam16_reconstruct(j-3) = temp(1);
        qam16_reconstruct(j-2) = temp(2);
        qam16_reconstruct(j-1) = temp(3);
        qam16_reconstruct(j) = temp(4);

    end


    for i = 1:8:length(qam16_reconstruct)-8
        qam16_reconstruct(i:i+7) = qam16_reconstruct(i+7 : -1: i);
    end

    s = num2cell(reshape(qam16_reconstruct,8,[])',2);
    b = cellfun(@(qam16_reconstruct) bin2dec(strrep(num2str(qam16_reconstruct),' ','')), s);
    img_qam=reshape(b,m,n);

    figure(4);
    subplot(2,1,1);
    imshow(img_qam, [0 255]);
    hold on;
    title("16-QAM Retrieved Signal, SNR = " + high_snr);


    low_snr = 5;
    qam16 = awgn(qam16, low_snr); %Low SNR
    qam16_demod = qamdemod(qam16, 16);

    qam16_reconstruct = zeros(1,length(stream));

    for i = 1:length(qam16_demod)
        j = i*4;
        temp = de2bi(qam16_demod(i),4,'left-msb');

        qam16_reconstruct(j-3) = temp(1);
        qam16_reconstruct(j-2) = temp(2);
        qam16_reconstruct(j-1) = temp(3);
        qam16_reconstruct(j) = temp(4);

    end


    for i = 1:8:length(qam16_reconstruct)-8
        qam16_reconstruct(i:i+7) = qam16_reconstruct(i+7 : -1: i);
    end

    s = num2cell(reshape(qam16_reconstruct,8,[])',2);
    b = cellfun(@(qam16_reconstruct) bin2dec(strrep(num2str(qam16_reconstruct),' ','')), s);
    img_qam=reshape(b,m,n);

    figure(4);
    subplot(2,1,2);
    imshow(img_qam, [0 255]);
    hold on;
    title("16-QAM Retrieved Signal, SNR = " + low_snr);
    
end

%Reference for converesion and reconstruction of image to bits and vice
%versa: https://stackoverflow.com/questions/47798711/how-to-transmit-image-to-bits-and-back-in-matlab
