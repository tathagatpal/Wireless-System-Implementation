image = imread('assign1.jpg');
stream = reshape((dec2bin(image,8)-'0').',1,[]);
[m, n] = size(image);

choice1 = input('Enter 1 for (7, 4) Hamming Code, 2 for Rate 1/2 Convolutional Code: ');

%(7, 4) HAMMING
if choice1 == 1
%(7, 4) Hamming Coding
stream_hamming = encode(stream,7,4);   % C(x) = A (7,4) Hamming codeword

choice = input('Enter 1 for QPSK, 2 for 64-QAM: ');
                 
%QPSK____________________________________________________
if choice == 1
    %QPSK Modulation WITH CHANNEL CODING
    qpsk_stream_encoding = zeros(1, length(stream_hamming)/2);
    for i = 1:2:length(stream_hamming)
        qpsk_stream_encoding(round(i/2)) = bin2dec(strcat(num2str(stream_hamming(i)),num2str(stream_hamming(i+1))));
    end
    qpsk_mod_coding = pskmod(qpsk_stream_encoding,4,pi/2);
    
    scatterplot(qpsk_mod_coding);
    hold on;
    title('Transmit Constellation with (7,4) Hamming Coding');
    
    %QPSK Modulation WITHOUT CHANNEL CODING 
    qpsk_stream = zeros(1, length(stream)/2);

    for i = 1:2:length(stream)
        qpsk_stream(round(i/2)) = bin2dec(strcat(num2str(stream(i)),num2str(stream(i+1))));
    end

    qpsk_mod = pskmod(qpsk_stream,4,pi/2);
    
    scatterplot(qpsk_mod);
    hold on;
    title('Transmit Constellation without Coding');
    figure;
    
    %RAYLEIGH CHANNEL WITH CHANNEL CODING__________________________________________________
    Eb_No = 0: 1: 30;
    snr_val = Eb_No + 10*log10(6);
    
    for ind = 1:1:length(snr_val)
        
        qpsk_mod_ray = qpsk_mod_coding; 

        %Rayleigh
        noise_var = 1/10^(snr_val(ind)/10);      %10.log(p/sigma^2)=SNR
        h = 1/2*randn(1,length(qpsk_mod_ray)) + 1j*randn(1,length(qpsk_mod_ray)); 
        % Send over Gaussian Link to the receiver
        y_out = h.*qpsk_mod_ray ...
            + sqrt(noise_var)*(randn(1,length(qpsk_mod_ray))+1j*randn(1,length(qpsk_mod_ray)));

        % (Ideal) Equalization to remove fading effects
        y_out = y_out./h;
        
        %QPSK Demodulation
        qpsk_demod = pskdemod(y_out, 4, pi/2);

        qpsk_reconstruct = zeros(1,length(stream));

        for i = 1:length(y_out)
            j = i*2;
            temp = de2bi(qpsk_demod(i),2,'left-msb');

            qpsk_reconstruct(j-1) = temp(1);
            qpsk_reconstruct(j) = temp(2);

        end

        %Hamming decoding
        hamming_decode = decode(qpsk_reconstruct,7,4);
        
        bit_error = 0;
        bit_total = 0;
        diff = stream - hamming_decode;
        bit_error = bit_error + sum(abs(diff));
        bit_total = bit_total + length(qpsk_mod_ray);
        BER_coding(ind) = bit_error / bit_total;

    end
    
    s = num2cell(reshape(hamming_decode,8,[])',2);
    b = cellfun(@(hamming_decode) bin2dec(strrep(num2str(hamming_decode),' ','')), s);
    img_qpsk=reshape(b,m,n);
    imshow(img_qpsk, [0 255]);
    figure;
    
    scatterplot(y_out);
    hold on;
    title('Received Constellation with (7,4) Hamming Coding');
    figure;
    
    
 %RAYLEIGH CHANNEL WITHOUT CHANNEL CODING_____________________________________________
    for ind = 1:1:length(snr_val)
        
        qpsk_mod_ray = qpsk_mod;

        %Rayleigh
        noise_var = 1/10^(snr_val(ind)/10);      %10.log(p/sigma^2)=SNR
        h = 1/2*randn(1,length(qpsk_mod_ray)) + 1j*randn(1,length(qpsk_mod_ray)); 

        % Send over Gaussian Link to the receiver
        y_out = h.*qpsk_mod_ray ...
            + sqrt(noise_var)*(randn(1,length(qpsk_mod_ray))+1j*randn(1,length(qpsk_mod_ray)));

        % (Ideal) Equalization to remove fading effects
        y_out = y_out./h;
        
        %QPSK Demodulation
        qpsk_demod = pskdemod(y_out, 4, pi/2);

        qpsk_reconstruct = zeros(1,length(stream));

        for i = 1:length(y_out)
            j = i*2;
            temp = de2bi(qpsk_demod(i),2,'left-msb');

            qpsk_reconstruct(j-1) = temp(1);
            qpsk_reconstruct(j) = temp(2);

        end

        %Hamming decoding        
        bit_error = 0;
        bit_total = 0;
        diff = stream - qpsk_reconstruct;
        bit_error = bit_error + sum(abs(diff));
        bit_total = bit_total + length(qpsk_mod_ray);
        BER_no_coding(ind) = bit_error / bit_total;
        
     end
    
    s = num2cell(reshape(qpsk_reconstruct,8,[])',2);
    b = cellfun(@(qpsk_reconstruct) bin2dec(strrep(num2str(qpsk_reconstruct),' ','')), s);
    img_qpsk=reshape(b,m,n);
    imshow(img_qpsk, [0 255]);
    
    scatterplot(y_out);
    hold on;
    title('Received Constellation without Coding');
    figure;
    
    %Theoretical Value
    EbN0Lin = 10.^(Eb_No/10);
    ber_theory = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+2)));
    
    semilogy(Eb_No, BER_coding...
        , 'r', 'linewidth', 1.5);
    hold on;
    semilogy(Eb_No, BER_no_coding...
        , 'b', 'linewidth', 1.5);
    semilogy(Eb_No, ber_theory...
        , 'k', 'linewidth', 1.5);
    xlabel('E_b/N_o (dB)'); ylabel('BER');
    title('BER v/s E_b/N_o (QPSK)');
    legend('With Convolutional Coding', 'Without Coding', 'Theoretical');
                 
                 
                 
%64-QAM____________________________________________________________________
  elseif choice == 2
    
    %WITH (7, 4) HAMMING CODING
    reversed_coding = stream_hamming(1:end-2);
    signal_qam_coding = [];
    for i = 1 : 6: length(reversed_coding) %taking 6 bits for modulation at a time
        signal_qam_coding = [signal_qam_coding 32*reversed_coding(i)+16*reversed_coding(i+1)+8*reversed_coding(i+2)+4*reversed_coding(i+3)+2*reversed_coding(i+4)+reversed_coding(i+5)];
    end

    qam_coding = qammod(signal_qam_coding, 64);
    
    scatterplot(qam_coding);
    hold on;
    title('Transmit constellation with (7,4) Hamming Code');
    
    
    %WITHOUT CODING
    reversed = stream(1:end-2);
    signal_qam = [];
    for i = 1 : 6: length(reversed) %taking 6 bits for modulation at a time
        signal_qam = [signal_qam 32*reversed(i)+16*reversed(i+1)+8*reversed(i+2)+4*reversed(i+3)+2*reversed(i+4)+reversed(i+5)];
    end

    qam_stream = qammod(signal_qam, 64);
    
    scatterplot(qam_stream);
    hold on;
    title('Transmit constellation without coding');
    figure;
   
    %RAYLEIGH CHANNEL WITH CHANNEL CODING__________________________________ 
    Eb_No = 5: 1: 30;
    snr_val = Eb_No + 10*log10(2);
    k = log2(64);
    
    for ind = 1:1:length(snr_val)
        
%         qam_coding_ray = qam_coding;
        qam_coding_ray = awgn(qam_coding, snr_val(ind),'measured');
        
        N = length(qam_coding_ray);
        %Rayleigh
        p = mean((abs(qam_coding_ray)).^2)/k;
        noise_var = sqrt(p/(2 * 10^(Eb_No(ind)/10)));
        noise_ray = noise_var*(randn(N,1) + 1j*randn(N,1));
        H = (1/sqrt(2)) * (randn(N,1)+1j *randn(N,1)); 
        y_out = H.*qam_coding_ray' + noise_ray;
        y_out = y_out./H;
        
        %QPSK Demodulation
        qam_demod=qamdemod(qam_coding_ray,64);
        qam_reconstruct=[];
        for i=qam_demod
            qam_reconstruct=[qam_reconstruct de2bi(i,6,'left-msb')];
        end

        qam_reconstruct=[qam_reconstruct 0 0];

        hamming_decode = decode(qam_reconstruct,7,4);
        
        bit_error = 0;
        bit_total = 0;
        diff = stream - hamming_decode;
        bit_error = bit_error + sum(abs(diff));
        bit_total = bit_total + length(qam_coding_ray);
        BER_coding(ind) = bit_error / bit_total;

    end
    
    s = num2cell(reshape(stream,8,[])',2);
    b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
    out=reshape(b,m,n);  
    imshow(out,[0,255]);
    hold on;
    title("Reconstructed using 64-QAM");
    
    scatterplot(y_out);
    title('Received Constellation with (7,4) Hamming Coding');
    figure;
    
    
    %RAYLEIGH CHANNEL WITHOUT CHANNEL CODING__________________________________ 
    
    for ind = 1:1:length(snr_val)
        
        qam_coding_ray = awgn(qam_stream, snr_val(ind),'measured');
                
        %Rayleigh
        p = mean((abs(qam_coding_ray)).^2)/k;
        noise_var = sqrt(p/(2 * 10^(Eb_No(ind)/10)));
        noise_ray = sqrt(noise_var)*(randn(1,length(qam_coding_ray))+1j*randn(1,length(qam_coding_ray)));
        H = 1/2*randn(1,length(qam_coding_ray)) + 1j*randn(1,length(qam_coding_ray)); 
        y_out = H.*qam_coding_ray + noise_ray;
        y_out = y_out./H;
        
        %QAM Demodulation
        qam_demod=qamdemod(qam_coding_ray,64);
        qam_reconstruct=[];
        for i = qam_demod
            qam_reconstruct=[qam_reconstruct de2bi(i,6,'left-msb')];
        end

        qam_reconstruct=[qam_reconstruct 0 0];
        
        bit_error = 0;
        bit_total = 0;
        diff = stream - qam_reconstruct;
        bit_error = bit_error + sum(abs(diff));
        bit_total = bit_total + length(qam_coding_ray);
        BER_no_coding(ind) = bit_error / bit_total;

    end
    
    s = num2cell(reshape(stream,8,[])',2);
    b = cellfun(@(stream) bin2dec(strrep(num2str(stream),' ','')), s);
    out=reshape(b,m,n);  
    imshow(out,[0,255]);
    hold on;
    title("Reconstructed using 64-QAM (without coding)");
    
    figure;
    scatterplot(y_out);
    title('Received Constellation without coding');
    figure;
    
    plot(Eb_No, BER_coding, 'r', 'linewidth', 1.5);
    hold on;
    title('BER vs E_b/N_o');
    plot(Eb_No, BER_no_coding, 'b', 'linewidth', 1.5);
    xlabel('E_b/N_o (dB)'); ylabel('BER');
%     ylim([0.003 1]); 
    legend('With (7, 4) Hamming Coding', 'Without Coding');
    
end
