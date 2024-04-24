%(7) Generalize program to have s1(t) and s2(t) as general waveforms or user defined.

%1-Simulation parameters
number_of_bits = 1e5;        
SNR_range_in_dB = 0:2:30;  
m = input('Enter the number of samples: ');   % number of samples to represent waveform
sampling_instant = input('Enter the number of sampling instant: ');             
s1 = input('Enter s1: ');     
s2 = input('Enter s2: ');   

%2-Generate random binary data vector
binary_data = randi([0 1],1,number_of_bits);           %Vector of  uniform random values, of size 1*e5

%3-Represent each bit with proper waveform, m = 20
waveForm = zeros(1,number_of_bits*m); %initialization of waveform, of size 20*e5
%we use the following for loop to Represent each bit with proper waveform
for i = 1:number_of_bits      
    if binary_data(i) == 1    
        waveForm((i-1)*m+1:i*m) = s1;     %waveform is filled with 20 position of ones 
    else         
        waveForm((i-1)*m+1:i*m) = s2;     %waveform is filled with 20 position of zeros  
    end
end

BER_simple = [];
BER_MF = [];
BER_Corr = [];

%Calculating threshold 
 V_Threshold = (sum(s1.*s1)-sum(s2.*s2))/2; 
 V_Threshold_simple = (s1+s2)/2;
    
for snr = 0:2:30
    
  %4-Add noise to samples
    noisyWF = awgn(waveForm, snr, 'measured');
    
 
  %5-comparing between simple detector, matched filter and correlator :
    
  %----------------------- simple detector --------------------------------
  
    simple_Received_TH = zeros(1,length(noisyWF)); %initialization
      
    for k = 1:length(noisyWF)
         if(noisyWF(k) >= V_Threshold_simple)
            simple_Received_TH(k) = 1;
        else
            simple_Received_TH(k) = 0;
        end 
    end
    
  %calculating simple detector BER
  [number, ratio] = biterr(waveForm, simple_Received_TH); %number, the number of bits that differ in the comparison, and ratio, the ratio of number to the total number of bits.
  BER_simple = [BER_simple ratio];
    
  %------------------------- matched filter -------------------------------
   
    h_mf = fliplr(s1-s2);      %reflection and shift with t=T "impulse response" of matched filter
    mf_output = zeros(1,number_of_bits); %initialization
    
    %we use the following for loop to performe convolution in bit-by-bit basis i.e. 20 samples by 20 samples
    for i = 1 : number_of_bits
        noisyWF_sample = noisyWF((i-1)*m+1:i*m);  %Extracting 20 samples
        mf_sample = conv(noisyWF_sample,h_mf);
        mf_output(i) = mf_sample(sampling_instant);  %convolution output at the sampling instant
    end

 received_TH = zeros(1,number_of_bits); %initialization
 
    %we use the following for loop to compare the output of MF with the theshould voltage 
    for j = 1:length(mf_output)
        if(mf_output(j) >= V_Threshold)
            received_TH(j) = 1;
        else
            received_TH(j) = 0;
        end 
    end
    
    %calculating matched filter BER
   [number, ratio] = biterr(binary_data, received_TH);
   BER_MF = [BER_MF ratio];
    
  %------------------- Correlator receiver -------------------------------- 
  
   g = s1 - s2;
   corr_output = zeros(1, number_of_bits); %initialization
    
  %we use the following for loop to operate the correlator fundamental (multiplication + integration /or summation)
    for i = 1:number_of_bits
        corr_output(i) = sum(noisyWF((i-1)*m+1:i*m).*g);      
    end
    
    corr_Received_TH = zeros(1,number_of_bits); %initialization
    
  %we use the following for loop to compare the output of correlator with the theshould voltage
    for q = 1:length(corr_output)
         if(corr_output(q) >= V_Threshold)
            corr_Received_TH(q) = 1;
        else
            corr_Received_TH(q) = 0;
        end 
    end
    
    %calculating correlator receiver BER
    [number, ratio] = biterr(binary_data, corr_Received_TH);
    BER_Corr = [BER_Corr ratio];
       
end

%plotting the BER curve against SNR 
figure;
semilogy(SNR_range_in_dB,BER_simple, SNR_range_in_dB, BER_MF, 'r',  SNR_range_in_dB, BER_Corr, 'o');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Simple detector', 'Matched Filter', 'Correlator');
title('Coherent Detection of Binary Signal');

% Calculate transmitted signal power
P_tx = mean(abs(waveForm).^2);
fprintf('Transmitted signal power : %f\n', P_tx);
