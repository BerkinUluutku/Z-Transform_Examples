function [alpha] = optimizedAlpha(signal)
%This function gives an 'optimized' alpha, ln(r), value for the modified
%Fourier transform. This function is written for increasing signals.

signal = signal - signal(1); 
signal = abs(signal); 
%Getting the absolute value because possible negative numbers mess up the 
%logarithm function. If any case a negative number appear, it should be due 
%to noise, which might as well changed to a positive number.



alpha = log(signal(2)/signal(end)) / ((2-length(signal)));


end

