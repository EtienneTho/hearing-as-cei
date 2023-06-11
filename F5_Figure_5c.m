% Etienne Thoret 2023 (c)
%
% Generate figure 5b: Separability of two pure tones - 2
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
% 
% This script is quite long to run, run only the first section and from
% line 83 if you just want to regenerate the figure from already processed
% data
%
% Figure_5c.csv - col 1: f0, col 2: number of sifting iterations

close all ;
clearvars;
addpath(genpath('./lib/'));

fs = 16000 ;
f_tab = (80:1000) ;

nFreq = length(f_tab) ;
nbAmp = 200 ;

theoretical_curve = zeros(1,nbAmp) ;
nbIMFs = 6 ;
nSiftingTab = (1:1000) ;

duration = .1 ;
t = linspace(0, duration, floor(duration*fs)) ;
a_tab = linspace(1,2,nbAmp) ;
tabDist = zeros(1,length(nSiftingTab)) ;
tabTheoretical = zeros(1,length(nSiftingTab)) ;
tabNmaxSifting = zeros(1,nFreq) ;
lgtSiftTab = length(nSiftingTab) ;
%%
% compute for one freq
i_f = 500;
f_tab(i_f)
for iSifting = 1:lgtSiftTab
    iSifting
    d_tab             = zeros(1,nbAmp);
    theoretical_curve = zeros(1,nbAmp);
    tt = (floor(.05*duration*fs):floor(.95*duration*fs)) ;
    parfor i_phi = 1:nbAmp
        theoretical_curve(i_phi) = plomp_(f_tab(i_f),f_tab(i_f)*a_tab(i_phi)) ;
        rng(randi(1000))
        phi1 = rand*2*pi ;
        %phi2 = rand*2*pi ;
        phi2 = 0 ;
        s1 = cos(2*pi*f_tab(i_f)*t) ;
        s2 = cos(2*pi*f_tab(i_f)*a_tab(i_phi)*t+phi2) ;
        x = (s1 + s2) ;

        imf = emdc_fix([],x, nSiftingTab(iSifting), nbIMFs) ;

        d_2 = sqrt(sum((imf(2,tt)-s1(tt)).^2))/sqrt(sum(s1(tt).^2)) ;

        d_tab(i_phi) = d_2 ; 
    end

        [~,index_dtab] = max(rescale(d_tab)) ;
        index_dtab = find(rescale(d_tab) < .98) ;
    tabDist(1,iSifting) = a_tab(min(index_dtab)) ;
end

%% compute the theoretical roughness curves for the different frequencies
for i_f2 = 1:nFreq
    for i_phi = 1:nbAmp
        theoretical_curve(i_phi) = plomp_(f_tab(i_f2),f_tab(i_f2)*a_tab(i_phi)) ;
    end
    theoretical_curve = theoretical_curve / max(theoretical_curve);

    [~,index_thc]  = max(theoretical_curve) ;
    tabTheoretical(i_f2) = a_tab(index_thc) ;
end

%% find the N sifting to optimize

for i_f3 = 1:nFreq
    tempInd = find(tabDist(:) > tabTheoretical(i_f3)) ;
    tabNmaxSifting(i_f3) = nSiftingTab(max(tempInd)) ;
end
save('./data/siftingOptimizationFramework.mat') ;
%%
load('./data/siftingOptimizationFramework.mat') ;
loglog(f_tab,tabNmaxSifting);
xlabel('f_0 (Hz)')
ylabel('Number of Sifting Iteration')
grid on
grid minor
opt.FontSize = 12;
axis([0,500,0,1000])
opt.BoxDim = [5, 5]; %[width, height]
title(['Influence of sifting'])
setPlotProp(opt);
matrixtowrite = [f_tab', tabNmaxSifting'] ;
% Figure_5c.csv - col 1: f0, col 2: number of sifting iterations
writematrix(matrixtowrite, './out/csv/Figure_5c.csv') ;

%%
saveas(gca,'./out/eps/F5_Figure_5c.eps','epsc')
saveas(gca,'./out/csv/F5_Figure_5c.fig','fig')

%% functions
function d2 = bindist(a, b)
    d2 = a>=b;
end

function tab = randspace(l,u,nb,s)
    rng(s) ;
    tab = l + (u-l)*rand(1,nb) ;
end


function [unique_x,y_mean] = sortScatterWithUnique(x,y)
    unique_x = unique(x) ;
    y_mean = zeros(1,length(unique_x)) ;
    for iX = 1:length(unique_x)
        x__ = x==unique_x(iX);
        if length(x__)>10
            y_mean(iX) = median(y(x==unique_x(iX))) ;
        end
    end
end