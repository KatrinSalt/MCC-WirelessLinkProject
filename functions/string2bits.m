%NO_PFILE
function bits = string2bits(str)
% STRING2BITS Converts an ASCII string to a sequence of bits.
% Will generate a column vector of logical values, most-significant-bit
% first. Output will be of length 8*length(string).
% Usage example:
% disp(string2bits('I <3 signal processing!');

% Convert ascii string to binary string
bits = reshape(dec2bin(str, 8).',1,[]);

% Convert to binary array
bits = bits - '0';

% Convert to logical data type
bits = logical(bits);
