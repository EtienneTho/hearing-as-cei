% Etienne Thoret 2023 (c)
%
% Generate figure 5b: Separability of two pure tones - 1
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%
% Figure_5b.csv - col 1: 1it, col 2: 10its, col 3: 100its, col 4: 1000its, col 5: theoretical roughness

close all ;
clearvars;
addpath(genpath('./lib/'));

fs = 16000 ;
f_tab = ones(1,1)*500 ;

nFreq = length(f_tab) ;
nbTrial = 10 ;
nbAmp = 100 ;

theoretical_curve = zeros(1,nbAmp) ;

duration = 1 ;
t = linspace(0, duration, floor(duration*fs)) ;
iterationTab = [1 10 100 1000] ;
d_tab = zeros(length(iterationTab),nbAmp) ;
a_tab = linspace(1,3,nbAmp) ;
densityTabTab = [] ;

nbIMFs = 6 ;
nbSiftingIterations = 1 ;

for i_f = 1:length(iterationTab) % loop on frequencies
    parfor_progress(nbAmp); % Initialize progress bar
    nbSiftingIterations = iterationTab(i_f) ;
    parfor i_phi = 1:nbAmp
        theoretical_curve(i_phi) = plomp_(f_tab,f_tab*a_tab(i_phi)) ;
        rng(randi(1000))
        phi1 = rand*2*pi ;
        phi2 = rand*2*pi ;
        phi2 = 0 ;

        s1 = cos(2*pi*f_tab*t) ;
        s2 = cos(2*pi*f_tab*a_tab(i_phi)*t+phi2) ;
        prod_s = 2*cos(pi*f_tab*(1+a_tab(i_phi))*t+phi2/2).*cos(pi*f_tab*(1-a_tab(i_phi))*t-phi2/2);
        x = (s1 + s2) ;

        y1 = s1;
        y2 = s2;

        imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;

        tt = (floor(.05*duration*fs):floor(.95*duration*fs)) ;
        d_1 = sqrt(sum((imf(1,tt)-prod_s(tt)).^2))/sqrt(sum(prod_s(tt).^2));
        d_1 = sqrt(sum((imf(1,tt)-x(tt)).^2))/sqrt(sum(x(tt).^2)/2) ;

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

for iTab = 1:length(iterationTab)
    plot(a_tab, rescale(d_tab(iTab,:)),'--');
    matrixtowrite = [matrixtowrite, rescale(d_tab(iTab,:))'];    
    hold on;
end

opt.FontSize = 12;
hold on;
theoretical_curve = (theoretical_curve(1,:)) ;
theoretical_curve = theoretical_curve / max(theoretical_curve) ;
plot(a_tab,theoretical_curve,':','Color','b')
legend('1 it.','10 its.','100 its.','1000 its.','Theoretical roughness')
xlabel('Interval')
ylabel('Dissonance (a.u.)')
opt.XLabel = 'Interval'; % xlabel
opt.YLabel = 'Separability (a.u.)'; %ylabel
grid on
grid minor
opt.BoxDim = [5, 5]; %[width, height]
title('Influence of the number of sifting iterations')
setPlotProp(opt);
matrixtowrite = [] ;
% Figure_5b_1.csv - col 1: 1it, col 2: 10its, col 3: 100its, col 4: 1000its, col 5: theoretical roughness
writematrix(matrixtowrite, './out/csv/Figure_5b.csv') ;


%%
saveas(gca,'./out/eps/F5_Figure_5b.eps','epsc')
saveas(gca,'./out/fig/F5_Figure_5b.fig','fig')

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