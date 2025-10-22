function [y] = sym_log(x)
y = sign(x).*log10(abs(x));
y(abs(x) < 0.5) = 0; 
end