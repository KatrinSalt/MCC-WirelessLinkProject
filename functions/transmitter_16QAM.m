function [tx,symbols] = transmitter_16QAM(bits_in,fs,Rb,preamble,const_var)

addpath('./functions');

if const_var == 1 % use qpsk instead
    const = [1+1i 1-1i -1-1i -1+1i]/sqrt(2);
else
    const = [1+3i 3+3i 3+1i 1+1i 1-1i 3-1i 3-3i 1-3i -3-1i -1-1i -1-3i -3-3i -3+3i -1+3i -1+1i -3+1i]/sqrt(2);
end
M=length(const);

bpsymb = ceil(log2(M));
Bsym = log2(M);  % bits/symb
Rsymb = Rb/Bsym; % Symbol rate (Rb = l*Rs, l = 2 for QPSK (2 bits per codeword))
Tsymb = 1/Rsymb; % Symbol time 
fsfd = floor(fs/Rsymb); % Number of samples per symbol 
% For rrc:
rolloff = 0.3;
span = 4;

%% TRANSMITTER

n = length(bits_in);
sym_num = n/bpsymb;

msg = buffer(bits_in,Bsym).';

% message -> symbol
symbols = const(bi2de(msg, 'left-msb')+1);


% add preamble and a delay 
pilot = 0.8*ones(1,300);
symbols_p = 0.7*[zeros(1,50) pilot 1.3*preamble symbols]; 


% upsamplig 
symbols_upsample = upsample(symbols_p, fsfd); % increases the sample rate of 'symbols' by inserting fsfd (1 bit) between samples => enable pulse shap. using conv.

% pulse convolution 
[pulse, ~] = rtrcpuls(rolloff, Tsymb, fs, span); 
tx = conv(pulse, symbols_upsample); 
tx = tx/max(abs(tx));
% figure, plot(real(tx))
% figure, pwelch(tx,[],[],[],fs,'centered','power'),title('TX signal transmitted')

end