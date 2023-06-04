% Etienne Thoret 2023 (c)
%
% Generate figure 4: Masking thresholds of pure tones in noise. 
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%

clearvars;

f_tab = 100 ;
Nf = length(f_tab) ;
fs = 16000 ;

nbF0 = 8 ;
freqTab = linspace(100,7500,nbF0);
freqTab = [250 500 1000 2000];

%%
nbF0 = length(freqTab) ;
nbBdw = 15 ;

nbLevels = 30 ;
nbTrials = 400 ;
nbReversals = 30;
nbLoop = 10 ;
duration = .1 ;
N = floor(fs*duration) ;
t = linspace(0,duration,N) ;

tab_inflexion = [] ;
levels = 10.^(linspace(10,30,nbLevels)/10) ;
maxIndTabTab = zeros(nbLoop,length(freqTab),nbBdw);
parfor_progress(nbLoop); % Initialize progress bar
for iLoop = 1:nbLoop
    iLoop
    maxIndTab = [];
    for freq = freqTab
        freq

        f_tab = freq ;
        bdw = linspace(1,2*freq,nbBdw) ;
        dist_sig_tab = zeros(1,1) ;
        dist_noise_tab = zeros(1,1) ;

        dist_sig_tab_repeat   = zeros(nbTrials,1) ;
        dist_noise_tab_repeat = zeros(nbTrials,1) ;
        tab_inflexion_tab = [] ;
    %     for iLevels = 1:nbLevels
            tab_inflexion = zeros(1,nbBdw) ;
            parfor iBdw = 1:nbBdw
                df = bdw(iBdw) ;
                bdw_min = max(f_tab-df/2,20) ;
                bdw_max = min(f_tab+df/2,fs/2-10) ;

                levelTest = 0 ;
                levelsTab = [] ;
                step = .3;
                detection_series = 0 ;
                iTrial = 1 ;
                reversals = [] ;
                while iTrial <= nbTrials
                    rng('shuffle')
                    noise = randn(1,N) ; %rand(1,N)*2-1; %
                    rng('shuffle')
                    phi   = 2*pi*rand ;
                    s1    = cos(2*pi*f_tab*t+phi)  ;
                    noise = fft_rect_filt(noise, bdw_min, min(bdw_max,fs/2-1), fs, 1) ;
                    %noise = cos(2*pi*min(bdw_max,fs/2-1)*t+phi)  ;   
                    
                    dist_sig    = emdMasker(s1, levelTest, noise, 0) ;
                    dist_notsig = emdMasker(s1, levelTest, noise, 1) ;
                    
                    dist_sig_tab = dist_sig ;
                    dist_noise_tab = dist_notsig ;

                    if iTrial > 1
                        if (dist_sig_tab < dist_noise_tab) ~= response
                           toSave = sum(s1.^2)/length(s1)/((sum((noise*levelTest).^2)/length(noise))) ;
                           reversals = [reversals toSave] ;
                        end
                    end
                    response = dist_sig_tab < dist_noise_tab ;

                    if response == 1 && detection_series == 3 
                        levelTest = levelTest + 3*step ;
                        detection_series = 0 ;
                    elseif response == 1 && detection_series ~= 3
                        detection_series = detection_series + 1 ;
                    else
                        levelTest = max(levelTest-step,0) ;
                        detection_series = 0 ;
                    end
                    if length(reversals) == nbReversals
                        break
                    end
                    levelsTab = [levelsTab levelTest] ;
                    iTrial = iTrial + 1 ;
%                     plot(levelsTab);drawnow
                end
                  tab_inflexion(iBdw) = nanmean(reversals) ;   
%                 df;
%                 tab_inflexion(1:iBdw);

            end
            tab_inflexion_tab = [tab_inflexion_tab; tab_inflexion] ;
            maxIndTab = [maxIndTab; tab_inflexion_tab] ;
    end
    
    maxIndTabTab(iLoop,:,:) = maxIndTab ;
    parfor_progress;
end

%%
Nplot = 1000 ;

matrixtowrite = linspace(0.1,2,Nplot)' ;
for iFreq = 1:nbF0
    iFreq
    x = linspace(0,2,2*nbBdw) ;
    y_mean = 20*log10(squeeze(maxIndTabTab(:,iFreq,:))) ;
    y_mean(y_mean==Inf) = NaN ;
    y_mean = nanmean(y_mean,1) ;
    
    y = [fliplr(y_mean) (y_mean)] ;
    p = polyfit(x,y,10);
    y1 = polyval(p,linspace(0,2,Nplot));
    semilogx(linspace(0.1,2,Nplot),y1-max(y1)) ;
    hold on;
    matrixtowrite = [matrixtowrite, (y1-max(y1))'];
end

writematrix(matrixtowrite, 'Figure_4.csv') ;

plot(linspace(0.1,2,Nplot),ones(1,Nplot)*-3,'--','Color','k')

xlabel('Interval')
ylabel('SNR for 80% correct (dB)')
legend('250 Hz', '500 Hz','1000 Hz', '2000 Hz','Location','northwest')
opt.XLabel = '\Delta f / f_0'; % xlabel

grid on
grid minor
axis([0,2,-15,1])
opt.BoxDim = [5, 5]; %[width, height]
title([''])
setPlotProp(opt);
%%
saveas(gca,'F4_Figure_4.eps','epsc2')
saveas(gca,'F4_Figure_4.fig','fig')

%%
function dist_ = emdMasker(sig, level, masker, notsig)        
        s1 = sig;
        
        if notsig == 0        
            x = s1 + masker * level  ;
        elseif notsig == 1
            x = masker ;
        end

        nbSiftingIterations = 1 ;
        nbIMFs = 6 ;         
        imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;            
        
        if size(imf,1) == 0
            imf(1,:) = 0;
        end
        
        % loop on imfs
        [sz,~] = size(imf);
        errDist = zeros(1,sz) ;
        for i = 1:sz
            errDist(i) = sqrt(sum((s1'-imf(i,:)').^2))/sqrt(sum(s1'.^2));
        end
        dist_ = min(errDist) ;
end
