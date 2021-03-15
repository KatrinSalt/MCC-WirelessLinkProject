function string = bits2string(bits)
 
string = char(bin2dec(reshape(char('0' + bits),8,[]).')).';

end