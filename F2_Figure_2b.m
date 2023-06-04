% Etienne Thoret 2023 (c)
%
% Generate figure 2b: Reconstruction of the missing fundamentals, averaged
% spectrum
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%

close all ;
clearvars;

%% generation of the harmonic complex
fs_initial = 16000 ; % sampling frequency
duration = 5 ; % duration of the sound
t = linspace(0, duration, floor(duration*fs_initial)) ; % times

signal = t ;
signal(:) = 0 ;

for freq = [1200 1600 2000]
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
opt.BoxDim = [5, 2]; %[width, height]

opt.YLim = [-100 0];
opt.XLim = [0 2200];
grid minor

matrixtowrite = [freqz, pow2db(sum(abs(pspectrum(IMFs_EMD_fix(tabIMF,1:end)',fs_initial)),2)), pow2db(abs(pspectrum(signal,fs_initial)))];
writematrix(matrixtowrite, 'Figure_2b.csv') ;

% apply
setPlotProp(opt);
saveas(gca,'F2_Figure_2b.eps','epsc')
saveas(gca,'F2_Figure_2b.fig','fig')


