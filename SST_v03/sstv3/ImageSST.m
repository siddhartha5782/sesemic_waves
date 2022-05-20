function ImageSST(fig)
%
% Display time-frequency result of the Synchrosqueezing transform ver 0.2
% Called by SSTguide.m
%
% INPUT:
%   fig.TFRtype = either CWT or SST
%   fig.TFR = sst.TFR;
%   fig.image = log(1+abs(srwt(1:40:end,:))); 
%   fig.t = sst.t(1:40:end);
%   fig.freq = freq;
%   fig.alpha = alpha;
%   fig.display = sst.display;
%
% OUTPUT:
%
% DEPENDENCY:
%   Wavelab (included), ContWavelet.m
%
% by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
% by Hau-tieng Wu 2012-03-03 (hauwu@math.princeton.edu)
%

if isfield(fig.display,'ytickadjust');
    fig.freq = fig.freq*fig.display.ytickadjust;
end

    %% plot time-frequency distribution, like CWT, STFT
if strcmp(fig.figtype,'TFD');

    [n, nscale] = size(fig.image);
    xtic = fig.t;

    if strcmp(fig.TFRtype, 'CWT');
    	ytic = linspace(2+(fig.TFR.oct-floor(log2(fig.TFR.scale))), log2(n)+2-floor(log2(fig.TFR.scale)), nscale);
    elseif strcmp(fig.TFRtype, 'STFT');
        ytic = fig.ytic;
    end

    imagesc(xtic, ytic, (fig.image)'); axis xy;

    if strcmp(fig.TFRtype, 'CWT') 
        ylabel('CWT. log_2(a)');
    elseif strcmp(fig.TFRtype, 'STFT');
        ylabel('STFT. \xi');
        %set(gca,'YTickLabel',fliplr(get(gca,'Ytick')));
    end
    %set(gca,'xtick',[]); 
    set(gca, 'fontsize', fig.display.fontsize);

end


    %% plot the synchrosqueezed TFD
if strcmp(fig.figtype,'SST')

    [n, nscale] = size(fig.image);
    fig.image   = fig.image(:,fliplr(1:1:nscale));
    fig.freq    = fig.freq(fliplr(1:1:nscale));

    xtic = fig.t;
    ytic = 1:nscale; 

	%% use quantile to reduce the possible extremal values to enhance the visualization
    QQ = quantile(fig.image(:), fig.display.quantile);
    fig.image(find(fig.image>QQ)) = QQ; clear QQ;

	%% enhance the visualization by locally adjust the contrast.
    if fig.display.SSTenhance
        for ii=1:n
            AX=fig.image(max(1,ii-20):min(size(fig.image,1),ii+20),:);
            amin = min(AX(:)); amax = max(AX(:));
            fig.image(ii,:) = ((fig.image(ii,:)-amin) ./ (amax-amin)) .*256;
        end
    end

    imagesc(xtic, ytic, fig.image'); axis('ij');

    if strcmp(fig.TFRtype, 'CWT') & fig.TFR.linear
        set(gca,'YTick',1:round(nscale/5):nscale);
    elseif strcmp(fig.TFRtype, 'CWT') & ~fig.TFR.linear
    	set(gca,'YTick',1:2*fig.TFR.nvoice:nscale);
    elseif strcmp(fig.TFRtype, 'STFT')
        set(gca,'Ytick',1:floor(nscale/5):nscale);
    end

    set(gca,'YTickLabel',round(fig.freq(get(gca,'YTick'))*100)/100);
    set(gca,'fontsize',fig.display.fontsize);
    ylabel('Synchrosqueezing. Hz');

end

    %% change the colormap
str = sprintf('colormap(1-%s(256))',fig.display.color);
eval(str)

end
