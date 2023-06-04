% Etienne Thoret 2023 (c)
%
% Generate figure 6: Dissonance
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%

clearvars;

nFreq = 12 ;% # of tested fequencies
f_tab = 2.^(linspace(0,12,nFreq)/12) * 261.63 ;


a_tab = [1 16/15 9/8 6/5 5/4 4/3 7/5 3/2 8/5 5/3 9/5 15/8 2] ;
a_tab = [a_tab linspace(1, 2.1, 100)];
a_tab = unique(sort(a_tab)) ;
nRatio = length(a_tab) ;
d = zeros(nFreq,nRatio) ; % initialize d => criteria on separation
R = zeros(nFreq,nRatio) ; % initialize d => criteria on separation
fs = 16000 ;

f_series = (1:15) ;
duration = .1 ;
t = linspace(0, duration, floor(duration*fs)) ;
nbRepeat = 50

parfor_progress(nFreq); % Initialize progress bar
for i_f = 1:nFreq % loop on frequencies
    for i_a = 1:nRatio % loop on amplitudes
        d_tab = zeros(1,nbRepeat) ; % local tab for the repetitions
        separation_tab = zeros(1,nbRepeat) ; % local tab for the repetitions

        parfor i_phi = 1:nbRepeat
            phi1 = rand*2*pi ;
            phi2 = rand*2*pi ;            

            s1 = sawtooth(2*pi*f_tab(i_f)*t+phi1) ;
            s2 = sawtooth(2*pi*f_tab(i_f)*a_tab(i_a)*t+phi2) ; 
            x = (s1 + s2) ;
            
            nbIMFs = 6 ;
            nbSiftingIterations = 1 ;
            nbSampleEnvelope = 1 ;

            imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;
            
            if size(imf,1) == 0
                imf(1,:) = 0;
                imf(2,:) = 0;
            elseif size(imf,1) == 1
                imf(2,:) = 0;
            end
            imf(isnan(imf)) = 0 ;

            imf1 = emdc([],s1);
            imf1 = emdc_fix([],s1, nbSiftingIterations, nbIMFs) ;
            if size(imf1,1) == 0
                imf1(1,:) = 0;
                imf1(2,:) = 0;
            elseif size(imf1,1) == 1
                imf1(2,:) = 0;
            end            
            imf1(isnan(imf1)) = 0 ;
            imf2 = emdc([],s2);
            imf2 = emdc_fix([],s2, nbSiftingIterations, nbIMFs) ;
            if size(imf2,1) == 0
                imf2(1,:) = 0;
                imf2(2,:) = 0;
            elseif size(imf2,1) == 1
                imf2(2,:) = 0;
            end
            imf2(isnan(imf2)) = 0 ;

            vecA =  sum(abs(fft(imf(:,:)'))')' ;            
            vecB = (sum(abs(fft(imf1(:,:)'))')' + sum(abs(fft(imf2(:,:)'))')') ;

            d_tab(i_phi) = mse(log10(vecA),log10(vecB))
        end
        d(i_f,i_a) =  nanmean(d_tab) ;
    end
    parfor_progress; 
end

parfor_progress(0); % Clean up
%%
d = 1-d;
%%
criterion = 10 ;
repeat = 1000 ;
densityTab = [];
for iF = 1:nFreq
    matTot = zeros(length(a_tab),length(a_tab)) ;
    for iRepeat = 1:repeat 
       mat_1 = tril(squareform(pdist(((d(iF,:))+randn(1,length(a_tab))/30)',@bindist))) ;
       mat_2 = tril(1 - mat_1)' ;
       matTot = matTot + mat_1 + mat_2;
    end
    density = btqn(matTot,criterion) ;
    densityTab = [density densityTab];
end

%% The Plomp & Levelt model

sweepingFrequency = 250:500;
sweepingFrequency = a_tab * 250 ;
consonanceCurve = zeros(length(sweepingFrequency),1);
spectrum = zeros(2,7);
for i = 1:length(sweepingFrequency)
    f = sweepingFrequency(i);
    for j = 1:7
        spectrum(1,j) = 250*(j+1);
        spectrum(2,j) = f*(j+1);
    end
    consonanceCurve(i) = plompSpectrum(spectrum);
end

plot(sweepingFrequency,1-consonanceCurve,'-b')
xlabel("Frequency")
ylabel("Dissonance")

%%
toPlot = mean(densityTab(:,:),2) ;
toPlot = (toPlot - min(toPlot)) ;
toPlot = toPlot / max(toPlot) ;

plot(a_tab, toPlot)
xticks([1 16/15 9/8 6/5 5/4 4/3 7/5 3/2 8/5 5/3 9/5 15/8 2])
xticklabels({'P1','m2', 'M2','m3', 'M3', 'P4','tt', 'P5', 'm6', 'M6', 'm7','M7', 'P8'})
set(gca, 'XScale', 'log')

opt.FontSize = 12;
hold on;
consonanceCurve = consonanceCurve / max(consonanceCurve) ;
plot(a_tab,consonanceCurve,'-.')

xlabel('Interval')
ylabel('Dissonance (a.u.)')
legend({'Simulated dissonance','Theoretical dissonance'})
opt.XLabel = 'Interval'; % xlabel
opt.YLabel = 'Dissonance (a.u.)'; %ylabel
grid on
grid minor
opt.BoxDim = [5, 2]; %[width, height]
setPlotProp(opt);
%%
matrixtowrite = [sweepingFrequency',toPlot,a_tab', consonanceCurve] ;
writematrix(matrixtowrite, 'Figure_6.csv') ;

%%
saveas(gca,'F6_Figure_6.eps','epsc')
saveas(gca,'F6_Figure_6.fig','fig')

%%


function d2 = bindist(a, b)
    d2 = a>=b;
end

%%

function plomp = plomp(f1,f2)
    fmin=min(f1,f2);
    fmax=max(f1,f2);
    s=0.24/(0.021*fmin+19.);
    plomp = (exp(-3.5*s*(fmax-fmin))-exp(-5.75*s*(fmax-fmin)));
end

function plompSpectrum = plompSpectrum(spectrum)
    [nSpectr, nFreq] = size(spectrum);

    c=0.0;
    for i = 1:nSpectr
        for j = (i+1):nSpectr
            for k = 1:nFreq
                for l = 1:nFreq
                    c=c+plomp(spectrum(i,k),spectrum(j,l));
                end
            end
        end
    end
    plompSpectrum = c;
end


