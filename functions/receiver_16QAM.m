function [bits_out, sym_rx] = receiver_16QAM(s_ch,fs,fc,Rb,preamble,bits_in,const_var)
%% RECEIVER: signal --> symbols --> message -->  bits 

addpath('./functions');

if const_var == 1 % use qpsk 
    const = [1+1i 1-1i -1-1i -1+1i]/sqrt(2);
else %16QAM
    const = [1+3i 3+3i 3+1i 1+1i 1-1i 3-1i 3-3i 1-3i -3-1i -1-1i -1-3i -3-3i -3+3i -1+3i -1+1i -3+1i]/sqrt(2);
end

M = length(const);
bpsymb = ceil(log2(M));
n = length(bits_in);
sym_num = n/bpsymb;

Bsym = log2(M);  % bits/symb
Rsymb = Rb/Bsym; % Symbol rate (Rb = l*Rs, l = 2 for QPSK (2 bits per codeword))
Tsymb = 1/Rsymb; % Symbol time 
fsfd = floor(fs/Rsymb); % Number of samples per symbol (is lower if bitrate is higher...)

% For rrc:
rolloff = 0.3;
span = 4;

%% MATCHED FILTER 
[MF, ~] = rtrcpuls(rolloff, Tsymb, fs, span); % matched filter 
MF_output = filter(MF, 1, s_ch); % run received signal through Matched Filter (MF)
MF_output = MF_output/max(abs(MF_output)); % normalize
%figure, pwelch(MF_output,[],[],[],fs,'centered','power')

%% CORRELATION PART 
pre_up = upsample(preamble,125); 
pre_up = 50*conv(pre_up, MF, 'same'); % 50 is only for seeing the difference better, does no ther change 
corr_new = conv(MF_output,fliplr(conj(pre_up)));
%figure, plot(abs(corr_new)), title('Correlation for frame syncronization')
% [corr_val, corr_pos] = max(abs(corrr_old));

%% frame sync
for cut = 0.5:0.5:1.5                   
     [y_val((cut-0.5)/0.5+1), ~] = max(abs(corr_new(cut*1e5:(cut+0.5)*1e5)));    % find maximum in each part of the correlation....     
end

temp_val = max(y_val) == abs(corr_new); % and choose the max of them
cor_pos_new = find(temp_val,1); % findind the x-position in correlation plot

MF_output_one = MF_output(cor_pos_new-2*(length(pre_up)+length(upsample(zeros(1,50),125))):cor_pos_new-2*(length(pre_up)+length(upsample(zeros(1,50),125)))+110000).'; % start form correlation minus premable etc and go to equal point + length of data 
%figure, plot(real(MF_output_one))

%% SYMBOL TIMING SYNCHRONIZATION 
% samples per symbol, find diffwerenc ebetween them and min should be at most open eye!
error_matrix = zeros(1, length(MF_output_one)-1);
[~, c] = size(error_matrix);
stepp = 3; % just how many we process, 3 equal every third data point. 
for i = 1:c-stepp
    idx1 = MF_output_one(:,i);
    idx2 = MF_output_one(:,i+stepp);
    error_matrix(i) = abs(idx2-idx1);
end 

% reshaped accoring to sampling size. creates vector with 125 columns for
% whole data and see where we should start sampling (1-125)
var = floor(length(error_matrix)/fsfd);
zer = fsfd - (length(error_matrix) - var * fsfd);
err_pad = [error_matrix, zeros(1, zer)];
err_reshaped = reshape(err_pad, [fsfd, length(err_pad)/fsfd]).';

[~, dem_point] = min(sum(err_reshaped)); % get demodulation point (aka best sample time after removed delay)
eyediagram(MF_output_one(dem_point:end), fsfd)
%% DOWNSAMPLING and fine frequency syncronization

MF_downsample_init = downsample(MF_output_one(dem_point:end), fsfd); %only initial value to use for sync!
%figure, plot(MF_downsample_init,'.'), title('Downsampled signal before fine freq. sync')

corr = conv(MF_downsample_init, fliplr(conj(preamble)));     % convolution with preamble pulse to find correlation 
%figure; plot(abs(corr), '.-r'),title('Correlation plot for preamble detection') % should have clear maximum
% find location of max corr
[~, Tmax] = max(abs(corr)); 
% plot(MF_downsample(Tmax+1:Tmax+nsymb)), '*')

% received PREAMBLE to use for freq sync
preamble_r = MF_downsample_init(Tmax - length(preamble)+1:Tmax);

% extra freq_sync
endpoint = 200;
for i=-endpoint:1:endpoint                                      % try to find a better value from initial place
    start =  i/1000;                                            % increase the offset little by little
    temp = preamble.*exp(1i*start*(0:length(preamble_r)-1));   	% calculate the signal for then
    corr = conv(preamble_r, fliplr(conj(temp)));                % calc correlation, aka when most accurate
    [cor_val(i+endpoint+1), ~] = max(abs(corr));                % find correlation 
    offset(i+endpoint+1) = start;                               % save offset 
end

[~,pos_] = max(cor_val); % find maximum correlation
figure, plot(offset,cor_val), title('Fine frequency syncronization'),xlabel('Frequency shift [Hz]'),ylabel('Correlation value') % to see so right point, should be a maximum here!
freq_offset = offset(pos_);

MF_downsample_final = MF_downsample_init.*exp(-1i*freq_offset*(0:length(MF_downsample_init)-1)); % shift back
figure, plot(MF_downsample_final,'.'), title('Downsampled signal after fine freq offset')

symbols_rx = MF_downsample_final(Tmax+1:Tmax+sym_num);
figure, plot(symbols_rx,'.'), title('Downsampled signal before phase syncronization') 
preamble_r_shifted = MF_downsample_final(Tmax - length(preamble)+1:Tmax);

pre_arr = [];
if const_var == 1 % QPSK
    for i=1:length(preamble)
        if real(preamble(i))>0 && imag(preamble(i))>0 % compare with one symbol (can choose any but best if preamble have many of that one.
            pre_arr = [pre_arr i];
        end
    end
else % 16 QAM
    for i=1:length(preamble)
        if (real(preamble(i))>0 && imag(preamble(i))>-1 && imag(preamble(i))<0) && real(preamble(i))<1 %%(real(preamble(i))<-1 && imag(preamble(i))<0 && imag(preamble(i))>-1)
            pre_arr = [pre_arr i];
        end
    end
end

phase_error = mean((angle(preamble_r_shifted(pre_arr)) - angle(preamble(pre_arr))).'); % find the phase error
user_data = symbols_rx.*exp(-1i*(phase_error));     % compensation of the phase error 

figure, plot(user_data,'.'), title('Downsampled signal after phase syncronization')

%% SCALING the signal
% if want other kind of scaling
% txEnergy = std(const); % Use the variation of the constellation 
% y_energy = std(user_data); % ...and the variation of the symbols...
% scale = txEnergy/y_energy; % find difference
scale = mean(abs(const))/mean(abs(user_data)); % scale with difference between input and output
sym_rx = user_data*scale;
figure, plot(const,'o'), hold on, plot(sym_rx,'.'), title('comparison sent and received symbol')

metric = abs(repmat(sym_rx.', 1, length(const)) - repmat(const, length(sym_rx), 1)).^2; % compute the distance to each possible symbol
% pick out the bits with shortest distances to the symbols and define to which quadrants the symbol belong to
[~, indx_r] = min(metric, [], 2); % find the closest for each received symbol
indx_r = indx_r.' - 1; 

% symbols --> message
msg_r = de2bi(indx_r, bpsymb, 'left-msb');

% message --> bits 
msg_r = msg_r.';
bits_out = msg_r(:).';
end