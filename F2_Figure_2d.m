% Etienne Thoret 2023 (c)
%
% Generate figure 2d: Combination tone, spectrum
% spectrum
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%
% Figure_2d.csv - col 1: CEI mean, col 2: CEI standard deviation

close all ;
clearvars;
addpath(genpath('./lib/'));

%% generation of the harmonic complex
fs_initial = 16000 ; % sampling frequency
duration = 5 ; % duration of the sound
t = linspace(0, duration, floor(duration*fs_initial)) ; % times

signal = t ;
signal(:) = 0 ;

for freq = [1000 1.4*1000]
    signal = signal + cos(2*pi*freq*t+rand*2*pi) ;
end

%% sifting 
nbIMFs = 6 ;
nbSiftingIterations = 1 ;
nbSampleEnvelope = 1 ;

IMFs_EMD_fix = emdc_fix([],signal, nbSiftingIterations, nbIMFs) ;

%% plots
figure
tabIMF = 1:nbIMFs;
[~,freqz]=pspectrum(signal,fs_initial);

plot(freqz,pow2db(sum(abs(pspectrum(IMFs_EMD_fix(tabIMF,1:end)',fs_initial)),2)));hold on;
plot(freqz,pow2db(abs(pspectrum(signal,fs_initial))));

addpath('./lib') ;
opt.XLabel = 'Frequency (Hz)'; % xlabel
opt.YLabel = 'Power (dB)'; %ylabel
opt.BoxDim = [5, 2.]; %[width, height]

opt.YLim = [-100 0];
opt.XLim = [500 1600];

% apply
setPlotProp(opt);


%%
saveas(gca,'./out/eps/F2_Figure_2d.eps','epsc')
saveas(gca,'./out/fig/F2_Figure_2d.fig','fig')

% Figure_2d.csv - col 1: CEI mean, col 2: CEI standard deviation
matrixtowrite = [freqz, pow2db(sum(abs(pspectrum(IMFs_EMD_fix(tabIMF,1:end)',fs_initial)),2)), pow2db(abs(pspectrum(signal,fs_initial)))];
writematrix(matrixtowrite, './out/csv/Figure_2d.csv') ;
