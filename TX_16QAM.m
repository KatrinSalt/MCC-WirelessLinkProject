%%% TRANSMITTER 
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
N = 10000;                              % Number of samples in a frame
frame_time = N/fs;                      % Time for 1 frame
time = (0:dt:dt*(N-1))';
RBW = 1/frame_time;
NFFT = 2^nextpow2(N);                   % Next power of 2 from length of y

fc = 20e6;                              % carrier frequency, iequal to fRF
Rb = 50e3;                              % bit rate (could have been changed)

const_var = 0;                          % zero means 16-QAM, 1 means QPSK ( working as flag)
%% message to send 
tx_str = ['Johnson traveled to Brussels on Wednesday for dinner with European Commission President Ursula von der Leyen. But the last-ditch effort failed to produce a breakthrough on thorny issues including fishing rights, government aid for companies and how disputes would be settled. '];
bits_in = string2bits(tx_str); % Convert the string to bits

% create preamble
const_QPSK = [1+1i 1-1i -1-1i -1+1i]/sqrt(2);
const_16QAM = [1+3i 3+3i 3+1i 1+1i 1-1i 3-1i 3-3i 1-3i -3-1i -1-1i -1-3i -3-3i -3+3i -1+3i -1+1i -3+1i]/sqrt(2);

if const_var == 1 % use qpsk instead
    preamble = const_QPSK(bi2de(get_preamble(const_var), 'left-msb')+1);
else
    preamble = const_16QAM(bi2de(get_preamble(const_var), 'left-msb')+1);
end

% transmitter function 
[send, symbols] = transmitter_16QAM(bits_in,fs,Rb,preamble,const_var);
s_tx = send.'; % this is the generated signal we send (real)
figure, pwelch(s_tx,[],[],[],fs,'centered','power'),xlim([-50 50]) % spectrum of transmitted signal

%% Setup the Tx
tx = comm.SDRuTransmitter(... 
'Platform','N200/N210/USRP2',...
'IPAddress','192.168.10.6',... % since that one is at our station
'CenterFrequency',20e6,...         % added freq noise here
'EnableBurstMode',1,...
'NumFramesInBurst',1,...
'InterpolationFactor',Interp_Factor,...
'MasterClockRate',MasterClock_Rate,...
'TransportDataType','int16');
%% Sending the signal  
currentTime = 0;
for k = 1:20000 % transmit many times...
    tx(s_tx);
    if k==1
        disp('transmitting')                % just to see that it has started
    end
    currentTime = currentTime+frame_time    % just to see that it is processing
end
release(tx);                                % must be don, otherwise cannot start transmission agian (for functioning of USRP)