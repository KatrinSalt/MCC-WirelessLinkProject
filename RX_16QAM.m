%%% RECEIVER 
%%% GROUP: ANTENNA AB 
clear all; close all; clc;
addpath('./functions');
%% This is the sampling rate for the digital mixer, do not change
MasterClock_Rate=100000000;
%% Interpolation factor for the Transmitter
Interp_Factor=64;
%% Decimation factor for the Receiver
Decimation_Factor=Interp_Factor;
%% Sampling rate and time
fs = MasterClock_Rate/Interp_Factor;    % sampling rate
dt = 1/fs;
N = 10000;
frame_time = N/fs;                      % Time for 1 frame
time = (0:dt:dt*(N-1))';
RBW = 1/frame_time;
NFFT = 2^nextpow2(N);                   % Next power of 2 from length of y

fc = 20e6;                              % carrier frequency, iequal to fRF
Rb = 50e3;                              % bit rate 

%% message to send 
tx_str = ['Johnson traveled to Brussels on Wednesday for dinner with European Commission President Ursula von der Leyen. But the last-ditch effort failed to produce a breakthrough on thorny issues including fishing rights, government aid for companies and how disputes would be settled. '];
bits_in = string2bits(tx_str); % Convert the string to bits

% create preamble
const_var = 0;      % 1= qpsk, 0=16-QAM
if const_var == 1
    const =  [1+1i 1-1i -1-1i -1+1i]/sqrt(2);
    preamble = const(bi2de(get_preamble(const_var), 'left-msb')+1);
else
    const = [1+3i 3+3i 3+1i 1+1i 1-1i 3-1i 3-3i 1-3i -3-1i -1-1i -1-3i -3-3i -3+3i -1+3i -1+1i -3+1i]/sqrt(2);
    preamble = const(bi2de(get_preamble(const_var), 'left-msb')+1);
end

%% Set up the RX
rx = comm.SDRuReceiver(...
    'Platform','N200/N210/USRP2',...
    'IPAddress','192.168.10.5',... % since that one is our station
    'CenterFrequency',20e6,... % center freq (could have used variable fc but we did not)
    'EnableBurstMode',1,...
    'NumFramesInBurst',1,...
    'DecimationFactor',Decimation_Factor,...
    'SamplesPerFrame',30*N,...
    'MasterClockRate',MasterClock_Rate,...
    'TransportDataType','int16');
%% Receiving the signal    
currentTime = 0;
m = 1; % how many times we run to get (but one was enough, if signal had been worse we could have sampled more versions)
for i=1:m
    [rx_data] = rx();
    rx_data = double(rx_data)/(2^8);            % scaling
    currentTime = currentTime+frame_time        % just for seeing that it is transmitting

    signal_rx = conj(rx_data.');                % conj back since HW has flipped to signal  
end
release(rx);                                    % must be done so the USRP can receive agian (finishes the operation)
%% plot input signal and its pectrum
%figure, plot(real(signal_rx)) % to plot received signal (good for seeing if we received correct things)
%figure, pwelch(signal_rx,[],[],[],fs,'centered','power'),xlim([-50
%50]),ylim([-100 40]), title('received signal') % to look at it powerspectrum (good for seeing if we received correct things)
%% CALCULATE COARSE FREQUENCY OFFSET
[pxx2, f2] = pwelch(signal_rx.',[],[],[],fs,'centered','power');
index2 = find(pxx2 == max(pxx2));
freq_shift_coarse = f2(index2)

receiver_input = signal_rx.*exp(2*-1i*pi*(freq_shift_coarse)*(0:length(signal_rx)-1)*1/fs); % shift back

% Plot spectrum of received and transmitted signal
figure(), subplot(2,1,1), pwelch(signal_rx,[],[],[],fs,'centered','power'),hold on, title('signal received')%,xlim([-50 50]),ylim([-100 50])
subplot(2,1,2),pwelch(receiver_input,[],[],[],fs,'centered','power'), title('signal after coarse shift back')%,xlim([-50 50]),ylim([-100 50]) % should be back around zero

% RECEIVER 
[bits_out,symbols_out] = receiver_16QAM(receiver_input.',fs,fc,Rb,preamble,bits_in,const_var); % our processing

% Received message
rx_str = bits2string(bits_out);

% Plotting constellation 
figure, plot(symbols_out,'*'), hold on, title('Symbols')
plot(const,'o'), axis equal, legend('RX','TX'), axis equal

% Calculating BER
bit_error = nnz(bits_out-bits_in); % should be 0  
% Comparison of transmitted and received message
disp(['Transmitted message:  ', tx_str])
disp(['Received message:     ', rx_str ]) 
disp(['Bit errors: ', num2str(bit_error)])