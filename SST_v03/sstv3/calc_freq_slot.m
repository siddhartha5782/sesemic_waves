function [idx] = calc_freq_slot(TFRtype, linear, iff, freq, alpha);

if strcmp(TFRtype, 'CWT');
    if linear
        idx = floor( (iff-freq(1))./alpha)+1;
    else
        idx = floor(log(iff/freq(1))./log(alpha))+1;
    end
else
    idx = floor((iff-freq(1))/alpha)+1;
end

