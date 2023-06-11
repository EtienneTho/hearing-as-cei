% Etienne Thoret 2023 (c)
%
% Generate figure 5a: Separability of two pure tones.
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%
% Figure_5a.csv - col 1:interval, col 2: 250 Hz, col 3: 500 Hz, col 4: 1000 Hz, col 5: 2000 Hz


close all ;
clearvars;
addpath(genpath('./lib/'));

fs = 16000 ;
addpath('./lib/') ;
f_tab = ones(1,1)*[250 500 1000 2000] ;

nFreq = length(f_tab) ;
nbTrial = 10 ;
nbAmp = 100 ;

theoretical_curve = zeros(nFreq,nbAmp) ;

duration = 1 ;
t = linspace(0, duration, floor(duration*fs)) ;
d_tab = zeros(nFreq,nbAmp) ;
a_tab = linspace(1,3,nbAmp) ;
densityTabTab = [] ;

nbIMFs = 6 ;
nbSiftingIterations = 1 ;

for i_f = 1:nFreq % loop on frequencies
    parfor_progress(nbAmp); % Initialize progress bar
    parfor i_phi = 1:nbAmp
        theoretical_curve(i_f,i_phi) = plomp_(f_tab(i_f),f_tab(i_f)*a_tab(i_phi)) ;
        rng(randi(1000))
        phi1 = rand*2*pi ;
        phi2 = rand*2*pi ;
        phi2 = 0 ;


        s1 = cos(2*pi*f_tab(i_f)*t) ;
        s2 = cos(2*pi*f_tab(i_f)*a_tab(i_phi)*t+phi2) ;
        prod_s = 2*cos(pi*f_tab(i_f)*(1+a_tab(i_phi))*t+phi2/2).*cos(pi*f_tab(i_f)*(1-a_tab(i_phi))*t-phi2/2);
        x = (s1 + s2) ;


        y1 = s1;
        y2 = s2;

        imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;

        tt = (floor(.05*duration*fs):floor(.95*duration*fs)) ;
        d_1 = sqrt(sum((imf(1,tt)-prod_s(tt)).^2))/sqrt(sum(prod_s(tt).^2));% ...        
        d_1 = sqrt(sum((imf(1,tt)-x(tt)).^2))/sqrt(sum(x(tt).^2)/2) ;%.* ...

        d_2 = abs(.5-sqrt(sum((imf(1,tt)-y2(tt)).^2))/sqrt(sum(y2(tt).^2))) ;
        d_2 = abs((sum((imf(2,tt)-y1(tt)).^2))/(sum(y1(tt).^2))) ;

        d_temp = d_1  ;

        d_tab(i_f,i_phi) = d_2 ; %(d_1 .* d_2) ;
        d_tab_2(i_f,i_phi) = d_2 ;
        parfor_progress;
        
    end
    d = mean(d_tab(:,:),1) ;
    d2 = mean(d_tab_2(:,:),1) ;

end

%%
matrixtowrite = a_tab' ;
for iTab = 1:nFreq
    plot(a_tab, rescale(d_tab(iTab,:)),'-');
    matrixtowrite = [matrixtowrite, rescale(d_tab(iTab,:))'];
    hold on;
end
% Figure_5a.csv - col 1:interval, col 2: 250 Hz, col 3: 500 Hz, col 4: 1000 Hz, col 5: 2000 Hz
writematrix(matrixtowrite, './out/csv/Figure_5a.csv') ;

opt.FontSize = 12;
legend('250 Hz', '500 Hz', '1000 Hz', '2000 Hz')
opt.XLabel = 'Interval'; % xlabel
opt.YLabel = 'Separability (a.u.)'; %ylabel
grid on
grid minor
opt.BoxDim = [5, 5]; %[width, height]
title('Influence of frequency')
setPlotProp(opt);

%%
saveas(gca,'./out/eps/F5_Figure_5a.fig','epsc')
saveas(gca,'./out/fig/F5_Figure_5a.fig','fig')

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