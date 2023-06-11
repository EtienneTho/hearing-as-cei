% Etienne Thoret 2023 (c)
%
% Generate figures 2a/b: Reconstruction of the missing fundamentals
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%

close all ;
clearvars;
addpath(genpath('./lib/'));

%% generation of the harmonic complex
filename = './data/voice_lowpass_filtered_300Hz.wav' ;
fs_target = 16000 ;
[signal,fs_initial] = audioread(filename) ;
signal = signal ./ 1.01 / max(abs(signal)) ;
signal = resample(signal,fs_target,fs_initial) ;
fs_initial = fs_target ;

%% sifting 
nbIMFs = 6 ;
nbSiftingIterations = 1 ;
nbSampleEnvelope = 1 ;

IMFs_EMD_fix = emdc_fix([],signal, nbSiftingIterations, nbIMFs) ;


%%
limSup = 2000 ; % in hz
wdw = 564*4 ;
hopsize = wdw/16 ;
[stEMD, stFourier] = plotEMDvsFourier(signal, wdw, hopsize, limSup, nbSiftingIterations, nbIMFs, fs_initial, fs_target);

%%
saveas(gca,'./out/eps/F2_Figure_2a_b.eps','epsc')
saveas(gca,'./out/fig/F2_Figure_2a_b.fig','fig')

%%
function [stEMD, st_fourier] = plotEMDvsFourier(signal, wdw,hopsize, limSup, nbSiftingIterations, nbIMFs, fs_initial, fs_target)
    
    x = signal ;
    nsegments = floor((length(x)-wdw) / hopsize) ;
    stEMD = [] ;
    st_fourier = [] ;
    for i = 1:nsegments
        fprintf('IMF %i %i\n',i,nsegments) ;
        segment = x((i-1)*hopsize+1:(i-1)*hopsize+1+wdw) ; % .* chebwin(wdw+1) ;
        IMFs_EMD_fix = emdc_fix([],segment, nbSiftingIterations, nbIMFs) ;
        pow_ = pow2db(sum(abs(pspectrum(IMFs_EMD_fix(:,1:end)',fs_initial)),2)) ;    
        stEMD = [stEMD (pow_)] ;
        pow_fourier = pow2db(abs(pspectrum(segment,fs_initial))) ;    
        st_fourier = [st_fourier (pow_fourier)] ;    
    end

    fig = figure(1);
    h = axes(fig,'visible','off'); 
    tiledlayout(2,1);
    
    nexttile
    colormap('parula')
    imagesc(flipud(stEMD(1:floor(limSup/fs_target*length(pow_)),:)),[-80 0]) ;
    title('Short-Term Cascaded Envelope Interpolation')    
    limSupFrequency = floor(limSup/fs_target*length(pow_)) ;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    duration = length(signal)/16000 ;
    frequency_sample = length(stEMD(:,1)) 
    xticks(floor([0 1 2 3] * length(stEMD(1,:)) / duration )+1)
    xticklabels({'0', '1', '2', '3'})
    yticks(floor([100 fs_target/4 fs_target/2 ] / (fs_target/2) * limSupFrequency) )
    yticklabels({fs_target/2, fs_target/4 0 })
    
    nexttile
    imagesc(flipud(st_fourier(1:floor(limSup/fs_target*length(pow_)),:)),[-80 0]) ;
    title('Short-Term Fourier Transform')    
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    xticks(floor([0 1 2 3] * length(stEMD(1,:)) / duration )+1)
    xticklabels({'0', '1', '2', '3'})    
    yticks(floor([100 fs_target/4 fs_target/2 ] / (fs_target/2) * limSupFrequency) )
    yticklabels({fs_target/2, fs_target/4 0 })    
    cb = colorbar;
    cb.Layout.Tile = 'east';
end


