%%% plot RECEIVER 
%%% GROUP: ANTENNA AB 
clear all; close all; clc;
%% This is the sampling rate for the digital mixer, do not change
MasterClock_Rate=100000000;
%% Interpolation factor for the Transmitter
Interp_Factor=64;
%% Decimation factor for the Receiver
Decimation_Factor=Interp_Factor;
%% Sampling rate and time
fs = MasterClock_Rate/Interp_Factor; %sampling rate
dt = 1/fs;
N = 10000;
frame_time = N/fs;% Time for 1 frame
time = (0:dt:dt*(N-1))';
RBW = 1/frame_time;
NFFT = 2^nextpow2(N); % Next power of 2 from length of y

fc = 20e6; % carrier frequency, iequal to fRF
Rb = 50e3; %b bit rate  (this is wahts hould have been changed...)

%% message to send 
tx_str = ['Johnson traveled to Brussels on Wednesday for dinner with European Commission President Ursula von der Leyen. But the last-ditch effort failed to produce a breakthrough on thorny issues including fishing rights, government aid for companies and how disputes would be settled. '];
bits_in = string2bits(tx_str); % Convert the string to bits

% create preamble
const_QPSK = [1+1i 1-1i -1-1i -1+1i]/sqrt(2);
const_16QAM = [1+3i 3+3i 3+1i 1+1i 1-1i 3-1i 3-3i 1-3i -3-1i -1-1i -1-3i -3-3i -3+3i -1+3i -1+1i -3+1i]/sqrt(2);

const_var = 0;          % qpsk =1, 16qam = 0
if const_var == 1 
    preamble = const_QPSK(bi2de(get_preamble(const_var), 'left-msb')+1);
else
    preamble = const_16QAM(bi2de(get_preamble(const_var), 'left-msb')+1);

end

%% Receiving the signal    
load ('100m_zero_2') % here is where to load the saved signal_rx data
%figure, plot(real(signal_rx))
% CALCULATE COARSE FREQUENCY OFFSET
[pxx2, f2] = pwelch(signal_rx.',[],[],[],fs,'centered','power');
index2 = find(pxx2 == max(pxx2));
freq_shift_coarse = f2(index2);
receiver_input = signal_rx.*exp(2*-1i*pi*(freq_shift_coarse)*(0:length(signal_rx)-1)*1/fs); % shift back
%% Plot spectrum of received and transmitted signal
figure(), subplot(2,1,1), pwelch(signal_rx,[],[],[],fs,'centered','power'),hold on, title('Received signal')%,xlim([-50 50]),ylim([-100 50])
subplot(2,1,2),pwelch(receiver_input,[],[],[],fs,'centered','power'), title('Received signal after coarse freq. shift back')%,xlim([-50 50]),ylim([-100 50]) % should be back around zero
%% RECEIVER 
%receiver_input = receiver_input/max(abs(receiver_input));
[bits_out,symbols_out] = receiver_16QAM(receiver_input.',fs,fc,Rb,preamble,bits_in,const_var);

% Received message
rx_str = bits2string(bits_out);

% Plotting constellation 
figure,
plot(symbols_out,'.'), hold on, title('Symbols received vs transmitted')
plot(const_16QAM,'o'), axis equal%, legend('RX','TX'),axis equal
%
% Calculating BER
bit_error = nnz(bits_out-bits_in); % should be  


%bit_error = abs((bits_out-ones(125,1)*bits_in));
%plot(sum(bit_error.'))
% Comparison of transmitted and received message
disp(['Transmitted message:  ', tx_str])
disp(['Received message:     ', rx_str ]) 
disp(['Bit errors: ', num2str(bit_error)])

