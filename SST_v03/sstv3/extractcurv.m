function [c] = extractcurv(curv);
%
% Curve extraction ver 0.2
%
% INPUT:
%   curv.TFRtype		= CWT or STFT
%   curv.stfd 		= the synchrosqueezed time-frequency domain
%   curv.freq 		= the discretzied freq range
%   curv.noctave 	= the total number of discretization in the frequency domain
%   curv.iff 		= the guessed instantaneous frequency
%   curv.range 		= the range around the guessed instantaneous frequency
%   curv.lambda 	= the penalty in curve extraction
%   curv.alpha 		= the step in the discretized freq range
%
% OUTPUT:
%   c: the curve related to the instantaneous frequency
%
% DEPENDENCY:
%   Wavelab (included), CurveExt.m
%
% by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
% by Hau-tieng Wu 2012-03-03 (hauwu@math.princeton.edu)
%


stfd = curv.stfd;

    %% get the frequency slot to speed up curve extraction
[idx] = calc_freq_slot(curv.TFRtype, curv.linear, curv.iff, curv.freq, curv.alpha);

%if strcmp(curv.TFRtype, 'CWT');
%    if curv.linear
%    	idx = floor( (curv.iff-curv.freq(1))./curv.alpha)+1;
%    else
%    	idx = floor(log(curv.iff/curv.freq(1))./log(curv.alpha))+1;
%    end
%else
%    idx = floor((curv.iff-curv.freq(1))/curv.alpha)+1;
%end

    %% to prevent overflow
tmp_stfd = abs(stfd(:,max(1,idx-curv.range):min(idx+curv.range,size(stfd,2))));

    %% call CurveExt_M.m to determine the parametrization of the dominant curve that can be visually seen in SST figure
c = CurveExt_M(tmp_stfd, curv.lambda);

    %% adjust the result
c = max(1, idx-curv.range)-1+c;

end
