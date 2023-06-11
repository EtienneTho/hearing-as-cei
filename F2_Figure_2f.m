% Etienne Thoret 2023 (c)
%
% Generate figure 2f: Tartini sounds
% spectrum
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%
% Figure_2f.csv - col 1: freq, col 2: cei, col 3: fourier

close all ;
clearvars;
addpath(genpath('./lib/'));

%% generation of the harmonic complex
fs = 16000 ;   % sampling frequency
duration = 3 ; % duration of the sound
t = linspace(0, duration, floor(duration*fs)) ; % times
rng(2)
freqzTab = (1:2000) ;
f0 = 400 ;
signal1 =  sawtooth_finite(2*pi*f0*t, 14) ;
signal2 =  sawtooth_finite(2*pi*f0*4/3*t, 14);

signal = signal1 + signal2;

fs_target = 16000 ;

fs = fs_target ;
duration = (length(signal)-1)/fs ;
start_ = 1;
end_ = 2;
duration = end_-start_
signal = signal(floor(start_*fs)+1:floor(end_*fs)) ;
signal = signal / 1.01 / max(abs(signal)) ;
t = linspace(0, duration, floor(duration*fs)) ;

% sifting
nbIMFs = 6 ;
nbSiftingIterations = 1 ;
nbSampleEnvelope = 1 ;

IMFs_EMD_fix = emdc_fix([],signal, nbSiftingIterations, nbIMFs) ;
IMFs_EMD = emdc([],signal) ;

tabIMF = 1:nbIMFs;

[~,freq] = pspectrum(IMFs_EMD(1,100:end-100)',fs);

pow_IMFs_EMD_fix = pow2db(sum(pspectrum(IMFs_EMD_fix(tabIMF,100:end-100)',fs),2)) ;

pow_Fourier = pow2db(pspectrum(signal)) ;

figure
freqzTab = (1:2000) ;
plot(freq(freqzTab),pow_IMFs_EMD_fix(freqzTab));hold on;
plot(freq(freqzTab),pow_Fourier(freqzTab));

addpath('./lib') ;
opt.XLabel = 'Frequency (Hz)'; % xlabel
opt.YLabel = 'Power (dB)'; %ylabel
opt.BoxDim = [5, 2]; %[width, height]
opt.YLim = [-90 0];
opt.XLim = [0 2000];
title('')
setPlotProp(opt);
%%
saveas(gca,'./out/eps/F2_Figure_2f.eps','epsc')
saveas(gca,'./out/fig/F2_Figure_2f.fig','fig')

% Figure_2f.csv - col 1: freq, col 2: cei, col 3: fourier
matrixtowrite = [freq(freqzTab), pow_IMFs_EMD_fix(freqzTab), pow_Fourier(freqzTab)];
writematrix(matrixtowrite, './out/csv/Figure_2f.csv') ;

function saw = sawtooth_finite(phi, N)
    saw = phi ;
    saw(:) = 0 ;
    for freq = (1:N)
        saw = saw + (-1)^(freq-1)*sin(phi*freq+rand*2*pi)/freq ;
    end
end
