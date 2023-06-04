% Etienne Thoret 2023 (c)
%
% Generate figure 3a: Cascaded envelope interpolation as an adaptive filter
% bank
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%

clearvars;

fs = 16000 ;

nbRepeat = 100 ;
duration = 1 ; 
N = floor(fs*duration) ;

nbIMFs = 6 ;
nbSiftingIterations = 1 ;
nbSampleEnvelope = 1 ;

imf_tab = zeros(nbRepeat,4096,nbIMFs) ;
cf_tab = zeros(1,nbRepeat*nbIMFs) ;
bdw_tab = zeros(1,nbRepeat*nbIMFs) ;
spec_x = [] ;

for i = 1:nbRepeat
    x = pinknoise(N)*40 ;    
    imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;            
   
    pspec = pspectrum(imf') ;
    pspec = pspec  ;
    
    for iIMFs = 1:nbIMFs
        [xc,lags] = xcorr(imf(iIMFs,:),x,1000) ;

        imf_tab(i,:,iIMFs) = pspectrum(xc)' ;
        imf_tab(i,:,iIMFs) = pspectrum(imf(iIMFs,:))' ;

        smoothed = smooth(pspec(:,iIMFs),1000) ;
        smoothed = smoothed / max(smoothed) ;
        [~,cf] = max(smoothed) ;
        [bdwMin,~] = find(smoothed(1:cf) < .5) ;
        bdwMin = max(bdwMin) ;
        [bdwMax,~] = find(smoothed(1:cf) > .5) ;
        bdwMax = max(bdwMax) ;
        bdw = bdwMax-bdwMin ;

        if isempty(bdw)
            bdw = 0 ;
        end
        bdw_tab((i-1)*nbIMFs+iIMFs) = bdw ;
        cf_tab((i-1)*nbIMFs+iIMFs) = cf ;
    end
    spec_x = [spec_x pspectrum(x(:))] ;
end


tabFilterMoy = zeros(nbIMFs,4096) ;
for iIMFs = 1:nbIMFs
    temp = [] ;
    for i = 1:nbRepeat
        temp = [temp ; imf_tab(i,:,iIMFs)] ;
    end    
    tabFilterMoy(iIMFs,:) = mean(temp,1);
end
%%
%subplot(141)
h = tiledlayout(1,4, 'TileSpacing', 'compact', 'Padding', 'none');
nexttile
spec = mean(spec_x,2) ;
[~,freq_] = pspectrum(x,fs) ;

plot(freq_,pow2db(spec),'-','Color','r','linewidth',.1) ;

hold on ;

plot(freq_,pow2db(tabFilterMoy'),'--','Color','k','linewidth',1) ;
set(gca, 'XScale', 'log') ;


set(gca,'xtick',[10 100 1000])
%set(gca,'xticklabel',[])
set(gca,'ytick',[-60 -40. -20. 0])
%set(gca,'yticklabel',[0 -20 -40])
axis('square')

axis([5 4000 -40 10]) ;
xlabel('Frequency (Hz)')
ylabel('Attenuation (dB)')
title('Pink noise')
legend('Signal spectrum','CEI filters','FontSize',6)

matrixtowrite = [freq_, pow2db(spec), pow2db(tabFilterMoy')];
writematrix(matrixtowrite, 'Figure_3a.csv') ;


%%
imf_tab = zeros(nbRepeat,4096,nbIMFs) ;
cf_tab = zeros(1,nbRepeat*nbIMFs) ;
bdw_tab = zeros(1,nbRepeat*nbIMFs) ;
spec_x = [] ;

for i = 1:nbRepeat
    i
    x =  sin(2*pi*100*linspace(0, duration, N))' + sin(2*pi*1000*linspace(0, duration, N))'  ;
    imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;            
    
%             imf = emdc([],x);    

    pspec = pspectrum(imf') ;
    pspec = pspec  ;
    
    for iIMFs = 1:nbIMFs
        [xc,lags] = xcorr(imf(iIMFs,:),x,1000) ;

        imf_tab(i,:,iIMFs) = pspectrum(xc)' ;
        imf_tab(i,:,iIMFs) = pspectrum(imf(iIMFs,:))' ;

        smoothed = smooth(pspec(:,iIMFs),1000) ;
        smoothed = smoothed / max(smoothed) ;
        [~,cf] = max(smoothed) ;
        [bdwMin,~] = find(smoothed(1:cf) < .5) ;
        bdwMin = max(bdwMin) ;
        [bdwMax,~] = find(smoothed(1:cf) > .5) ;
        bdwMax = max(bdwMax) ;
        bdw = bdwMax-bdwMin ;

        if isempty(bdw)
            bdw = 0 ;
        end
        bdw_tab((i-1)*nbIMFs+iIMFs) = bdw ;
        cf_tab((i-1)*nbIMFs+iIMFs) = cf ;
    end
    spec_x = [spec_x pspectrum(x(:))] ;
end


tabFilterMoy = zeros(nbIMFs,4096) ;
for iIMFs = 1:nbIMFs
    temp = [] ;
    for i = 1:nbRepeat
        temp = [temp ; imf_tab(i,:,iIMFs)] ;
    end    
    tabFilterMoy(iIMFs,:) = mean(temp,1);
end
%%
%subplot(142)
nexttile
spec = mean(spec_x,2) ;
[~,freq_] = pspectrum(x,fs) ;

plot(freq_,pow2db(spec),'-','Color','r','linewidth',.1) ;

hold on ;
plot(freq_,pow2db(tabFilterMoy'),'--','Color','k','linewidth',1) ;
set(gca, 'XScale', 'log') ;



set(gca,'xtick',[10 100 1000])
%set(gca,'xticklabel',[])
set(gca,'ytick',[-60 -40. -20. 0])
%set(gca,'yticklabel',[0 -20 -40])
axis('square')

axis([5 4000 -40 10]) ;
xlabel('Frequency (Hz)')
ylabel('Attenuation (dB)')
title('Pure tones')

matrixtowrite = [freq_, pow2db(spec), pow2db(tabFilterMoy')];
writematrix(matrixtowrite, 'Figure_3b.csv') ;

%%
imf_tab = zeros(nbRepeat,4096,nbIMFs) ;
cf_tab = zeros(1,nbRepeat*nbIMFs) ;
bdw_tab = zeros(1,nbRepeat*nbIMFs) ;
spec_x = [] ;

for i = 1:nbRepeat
    i
    x =  pinknoise(N)*40 + sin(2*pi*100*linspace(0, duration, N))' + sin(2*pi*1000*linspace(0, duration, N))'  ;
    imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;            
    
    pspec = pspectrum(imf') ;
    pspec = pspec  ;
    
    for iIMFs = 1:nbIMFs
        [xc,lags] = xcorr(imf(iIMFs,:),x,1000) ;

        imf_tab(i,:,iIMFs) = pspectrum(xc)' ;
        imf_tab(i,:,iIMFs) = pspectrum(imf(iIMFs,:))' ;

        smoothed = smooth(pspec(:,iIMFs),1000) ;
        smoothed = smoothed / max(smoothed) ;
        [~,cf] = max(smoothed) ;
        [bdwMin,~] = find(smoothed(1:cf) < .5) ;
        bdwMin = max(bdwMin) ;
        [bdwMax,~] = find(smoothed(1:cf) > .5) ;
        bdwMax = max(bdwMax) ;
        bdw = bdwMax-bdwMin ;

        if isempty(bdw)
            bdw = 0 ;
        end
        bdw_tab((i-1)*nbIMFs+iIMFs) = bdw ;
        cf_tab((i-1)*nbIMFs+iIMFs) = cf ;
    end
    spec_x = [spec_x pspectrum(x(:))] ;
end


tabFilterMoy = zeros(nbIMFs,4096) ;
for iIMFs = 1:nbIMFs
    temp = [] ;
    for i = 1:nbRepeat
        temp = [temp ; imf_tab(i,:,iIMFs)] ;
    end    
    tabFilterMoy(iIMFs,:) = mean(temp,1);
end
%%
%subplot(143)
nexttile
spec = mean(spec_x,2) ;
[~,freq_] = pspectrum(x,fs) ;

plot(freq_,pow2db(spec),'-','Color','r','linewidth',.1) ;

hold on ;
plot(freq_,pow2db(tabFilterMoy'),'--','Color','k','linewidth',1) ;
set(gca, 'XScale', 'log') ;

set(gca,'xtick',[10 100 1000])
set(gca,'ytick',[-60 -40. -20. 0])
axis('square')

axis([5 4000 -40 10]) ;
xlabel('Frequency (Hz)')
ylabel('Attenuation (dB)')
title('Pink noise + Pure tones')

matrixtowrite = [freq_, pow2db(spec), pow2db(tabFilterMoy')];
writematrix(matrixtowrite, 'Figure_3c.csv') ;

%%
load('./data/results_efficientCoding_all_speech.mat')
%subplot(144)
nexttile
cfTab_vec = (100:50:8000) ;
cfTab_vec = logspace(2,5,20); 

bdwTab_mean = zeros(1,length(cfTab_vec)) ;
bdwTab_std  = zeros(1,length(cfTab_vec)) ;
bdwTab_nb   = zeros(1,length(cfTab_vec)) ;

for i = 1:(length(cfTab_vec)-1)
    ind = find(cf_tab < cfTab_vec(i+1) & cf_tab > cfTab_vec(i)) ;
    bdwTab_mean(i) = mean(bdw_tab(ind)) ;
    bdwTab_std(i)  = std(bdw_tab(ind)) ;
    bdwTab_nb(i)   = length(bdw_tab(ind)) ;
end

vec = (0:8000) ;
ERB_linear = 24.7 * (4.37 *  vec / 1000 + 1) ;

ERB_poly = (6.23 * (vec/1000).^2 + 93.39 * (vec/1000) + 28.52) ;
Bark = 6 * asinh((vec)/600);

idxNan = isnan(log10(cf_tab));
cf_tab(idxNan) = [] ;
bdw_tab(idxNan) = [] ;
idxNan = isnan(log10(bdw_tab));
cf_tab(idxNan) = [] ;
bdw_tab(idxNan) = [] ;
idxNan = cf_tab < 100 ;
cf_tab(idxNan) = [] ;
bdw_tab(idxNan) = [] ;
idxNan = cf_tab > 7000 ;
cf_tab(idxNan) = [] ;
bdw_tab(idxNan) = [] ;

fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,-10]);
ft = fittype(@(a, b, x) a*x+b,'options',fo) ;

[fit_] = fit(log10(cf_tab),log10(bdw_tab),ft) ;

fitted_curve = fit_.a * log10(cfTab_vec(1:end-1))'+fit_.b  ;
scatter(cf_tab,bdw_tab,1.2, 'filled');hold on;
axis([0 3000 0 1000])
set(gca,'yscale','log')
set(gca,'xscale')
onethird_linear = .27 * cfTab_vec(1:end-1)' ;
ERB_linear = 24.7 * (4.37 *  cfTab_vec(1:end-1)' / 1000 + 1) ;
errorbar(cfTab_vec,bdwTab_mean,bdwTab_std,'lineWidth',1','color','black');hold on;
plot(cfTab_vec(1:end-1)', ERB_linear,'lineWidth',2);hold on;
plot(cfTab_vec(1:end-1)', onethird_linear,'lineWidth',2);hold on;
xlabel('Center Frequency (Hz)')
ylabel('Bandwidth (Hz)')
legend({'CEI','CEI_{average}','ERB','1/3 octave'},'Location','southeast','FontSize',8)
axis([100 4000 1 1500])
grid on
hold on;
axis square
title('Complex sounds')


matrixtowrite = [cf_tab, bdw_tab];
writematrix(matrixtowrite, 'Figure_3d_1.csv') ;

matrixtowrite = [cfTab_vec(1:end-1)', ERB_linear, onethird_linear];
writematrix(matrixtowrite, 'Figure_3d_2.csv') ;

matrixtowrite = [bdwTab_mean', bdwTab_std'];
writematrix(matrixtowrite, 'Figure_3d_3.csv') ;

%%
saveas(gca,'F3_Figure_3a.eps','epsc2')
saveas(gca,'F3_Figure_3a.fig','fig')

% matrixtowrite = [freq(freqzTab), pow_IMFs_EMD_fix(freqzTab), pow_Fourier(freqzTab)];
% writematrix(matrixtowrite, 'Figurre_2d.csv') ;
