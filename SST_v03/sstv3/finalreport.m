%
% final report of the Synchrosqueezing transform
% called by synchrosqueezing.m (2012-03-03)
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	%% display results
	%% detemine the time scale. By default, the input time unit is second,
	%% and second/minute/hour time scales are adjusted accordingly
	%% Or you can use your own time scale
if strcmp(sst.display.timescale,'hr')
    timescale=60*60;
elseif strcmp(sst.display.timescale,'min')
    timescale=60;
elseif strcmp(sst.display.timescale,'sec')
    timescale=1;
elseif strcmp(sst.display.timescale,'user');
    timescale=1;
    fprintf('Use the user-defined time scale.\n');
else
    timescale=1;
    fprintf('(DEBUG) No correct timescale is assigned. Use second.\n');
end

Xno = 1 + sst.display.dominant; 
Yno = sst.display.signal + sst.display.TFD + sst.display.SST + sst.display.component + sst.display.AM + sst.display.phase + sst.display.iff; 
setupfig0;

for qq = 1:Xno*Yno
    eval(['rslt.h',num2str(qq),'=h',num2str(qq),';']); 
end

Kno=1;

if sst.display.signal
    eval(['subplot(h',num2str(Kno),');']); 
    eval(['rslt.h',num2str(Kno),'=h',num2str(Kno),';']); Kno=Kno+1;
    plot(sst.t/timescale, x);
    set(gca,'fontsize',sst.display.fontsize); %set(gca,'xtick',[]); 
    axis tight;
    if Xno==2; 
	eval(['subplot(h',num2str(Kno),');']); axis off; Kno=Kno+1;
    end
end


if sst.display.TFD
    eval(['subplot(h',num2str(Kno),');']); 
    %eval(['rslt.h',num2str(ceil(Kno/Xno)),'=h',num2str(Kno),';']); 
    Kno = Kno+1;

    if sst.display.log
        fig.image = log(1+abs(rslt.tfd(1:sst.display.resolution:end,:))); 
    else
        fig.image = abs(rslt.tfd(1:sst.display.resolution:end,:)); 
    end

    fig.t 	= sst.t(1:sst.display.resolution:end)/timescale;
    fig.alpha 	= alpha;
    fig.TFRtype = sst.TFRtype;
    fig.TFR     = sst.TFR;
    fig.display = sst.display;
    fig.figtype = 'TFD';
    fig.ytic 	= tfd_ytic./sst.TFR.rescale;

    ImageSST(fig); clear fig;
  
    if Xno == 2; 
	eval(['subplot(h',num2str(Kno),');']); axis off; Kno=Kno+1;
    end

end

if sst.display.SST
    eval(['subplot(h',num2str(Kno),');']); 
    %eval(['rslt.h',num2str(ceil(Kno/Xno)),'=h',num2str(Kno),';']); 
    Kno=Kno+1;

    if sst.display.log
        fig.image = log(1+abs(rslt.stfd(1:sst.display.resolution:end,:))); 
    else
        fig.image = abs(rslt.stfd(1:sst.display.resolution:end,:)); 
    end

    fig.t 	= sst.t(1:sst.display.resolution:end)/timescale;
    fig.freq 	= freq./sst.TFR.rescale;
    fig.alpha 	= alpha;
    fig.TFRtype = sst.TFRtype;
    fig.TFR     = sst.TFR;
    fig.display = sst.display;
    fig.figtype = 'SST';

    ImageSST(fig); clear fig; 

	%% plot the extracted curve
    if sst.display.SSTcurv & sst.reconstruction
	hold on;
        for mm = 1:sst.extractcurv.no; 
            eval(['plot(sst.t/timescale, length(freq)-rslt.c',num2str(mm),'+1, ''g--'',''linewidth'',0.8);']);
        end; 
    end

	%% plot the band around the extracted curve (the band for reconstruction)
    if sst.display.reconband & sst.reconstruction
	hold on;
        for qq = 1 : sst.extractcurv.no

	    calc_band;

	    if qq == 1;
                eval(['plot(sst.t/timescale, length(freq)-rslt.c',num2str(qq),'+Z1, ''r--'',''linewidth'',0.5);']);
                eval(['plot(sst.t/timescale, length(freq)-rslt.c',num2str(qq),'-Z2, ''r--'',''linewidth'',0.5);']);
            else
                eval(['plot(sst.t/timescale, length(freq)-rslt.c',num2str(qq),'+Z1, ''g--'',''linewidth'',0.5);']);
                eval(['plot(sst.t/timescale, length(freq)-rslt.c',num2str(qq),'-Z2, ''g--'',''linewidth'',0.5);']);
            end

	end
    end


	%% plot the band around the guessed iff (the band for curve extraction)
    if sst.display.curvband & sst.reconstruction
        hold on;
        for qq = 1 : sst.extractcurv.no
            eval(['[idx] = calc_freq_slot(sst.TFRtype, sst.TFR.linear, sst.extractcurv.iff',num2str(qq),', rslt.freq, rslt.alpha);']);
            eval(['plot(sst.t/timescale, length(freq)-idx+',num2str(sst.extractcurv.range),', ''c--'',''linewidth'',0.5);']);
            eval(['plot(sst.t/timescale, length(freq)-idx-',num2str(sst.extractcurv.range),', ''c--'',''linewidth'',0.5);']);
        end
    end


	%% plot the "dominant frequency"
    if Xno == 2; 
	ytic = linspace(2+(sst.TFR.oct-floor(log2(sst.TFR.scale))),log2(length(sst.t))+2-floor(log2(sst.TFR.scale)),length(freq));
	eval(['subplot(h',num2str(Kno),');']);
	plot( abs(sum( rslt.stfd , 1 )), ytic); 
	set(gca,'ytick',[]);
	Kno = Kno+1;
    end; 

end


if sst.display.component
    eval(['subplot(h',num2str(Kno),');']);
    %eval(['rslt.h',num2str(Kno),'=h',num2str(Kno),';']); 
    Kno=Kno+1;

	%% plot the reconstructed components
    for qq = 1:sst.extractcurv.no
	hold on;
        eval(['plot(sst.t/timescale, real(rslt.recon',num2str(qq),')+2.5*(qq-1)*max(real(rslt.recon',num2str(qq),')),''r'');']);
    end

    ylabel('component');
    axis tight;

    if Xno==2; 
	eval(['subplot(h',num2str(Kno),');']); axis off; Kno = Kno+1; 
    end
end


if sst.display.AM
    eval(['subplot(h',num2str(Kno),');']);
    %eval(['rslt.h',num2str(Kno),'=h',num2str(Kno),';']); 
    Kno=Kno+1;

	%% plot the estimated AM
    for qq = 1:sst.extractcurv.no
        hold on;
        eval(['plot(sst.t/timescale, abs(rslt.recon',num2str(qq),')+2.5*(qq-1)*max(abs(rslt.recon',num2str(qq),')),''r'');']);
    end

    ylabel('AM');
    axis tight;

    if Xno==2; 
	eval(['subplot(h',num2str(Kno),');']); axis off; Kno=Kno+1;
    end

end


if sst.display.phase
    eval(['subplot(h',num2str(Kno),');']);
    %eval(['rslt.h',num2str(Kno),'=h',num2str(Kno),';']); 
    Kno=Kno+1;

        %% plot the estimated phase
    for qq = 1:sst.extractcurv.no
        hold on;
        eval(['plot(sst.t/timescale, phase(rslt.recon',num2str(qq),'./abs(rslt.recon',num2str(qq),')),''r'');']);
    end

    ylabel('phase');
    axis tight;

    if Xno==2; 
	eval(['subplot(h',num2str(Kno),');']); axis off; Kno = Kno+1; 
    end
end

if sst.display.iff
    eval(['subplot(h',num2str(Kno),');']);
    %eval(['rslt.h',num2str(Kno),'=h',num2str(Kno),';']); 
    Kno=Kno+1;

        %% plot the estimated phase
    for qq = 1:sst.extractcurv.no
        hold on;
        eval(['plot(sst.t/timescale, freq(rslt.c',num2str(qq),'),''r'');']);
    end

    ylabel('iff');
    axis tight;

    if Xno==2; 
	eval(['subplot(h',num2str(Kno),');']); axis off; Kno=Kno+1;
    end
end


eval(['subplot(h',num2str(Kno-Xno),');']);
set(gca,'xtickMode', 'auto');
xlabel(['time (',sst.display.timescale,')']);

if strcmp(sst.display.timescale,'user');

    set(gca, 'xtick', sst.display.xtick);
    set(gca, 'xticklabel', sst.display.xticklabel);
    xlabel(['time (',sst.display.xlabel,')']);

end

if sst.display.mesh

    figure;
    mesh(rslt.freq, sst.t(100:2000), abs(stfd(100:2000,:)));
    axis tight;
    str = sprintf('colormap(1-%s(256))','hot');
    eval(str)
    ylabel('time'); xlabel('freq');

end

