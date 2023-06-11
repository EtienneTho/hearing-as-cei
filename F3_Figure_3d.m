% Etienne Thoret 2023 (c)
%
% Generate figure 3d: Cascaded envelope interpolation as an adaptive filter
% bank
%
% It computes and saves the results for the panel d in './data/results_efficientCoding_all_speech.mat'
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%

close all ;
clearvars;
addpath(genpath('./lib/'));

nbRepeat = 1500 ;

nbIMFs = 6 ;
nbSiftingIterations = 1 ;

imf_tab = zeros(nbRepeat,4096,nbIMFs) ;

cf_tab = [] ;
bdw_tab = [] ;
pow_tab = [] ;

path_ = './cmos/all/' ;
addpath(genpath(path_)) ;
filelist = dir([path_ '*.wav']) ;
nbChunk = 0 ;
fs_target = 16000 ;
for i = 1:nbRepeat
    i
    [sig,fs] = audioread(filelist(i).name) ;
    sig = resample(sig,fs_target,fs) ;
    durationChunck = 800 ;
    nbChuncks = length(sig) / durationChunck  ;
    for iChunk = 1:nbChuncks-1
       nbChunk = nbChunk + 1;
       x = sig((iChunk-1)*durationChunck + 1:(iChunk)*durationChunck + 1) ;
       imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;            

        imf(isnan(imf)) = 0 ;

        [pspec,freqz] = pspectrum(imf',fs_target) ;

        nbIMFs_ = size(imf);
        nbIMFs_ = nbIMFs_(1) ;
        for iIMFs = 1:min(nbIMFs,nbIMFs_)

            [xc,lags] = xcorr(imf(iIMFs,:),x,1000,'normalized') ;
            smoothed = pspec(:,iIMFs) ;
            pow = sum(smoothed) ;
            pow_tab = cat(1,pow_tab, pow) ;

            smoothed = smoothed / max(smoothed) ;                        
            cf = sum(smoothed .* freqz) / sum(smoothed) ;
            bdw = sqrt(1/length(smoothed) * sum(smoothed.*(freqz-cf).^2)) ;

            if isempty(bdw)
                bdw = 0 ;
            end

            bdw_tab = cat(1,bdw_tab, bdw) ;
            cf_tab = cat(1,cf_tab, cf) ;
        end
    end

end



%%
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
legend({'CEI','CEI_{average}','ERB','1/3 octave'},'Location','southeast')
axis([100 4000 1 1500])
grid on
hold on;
axis square

%saveas(gca,'F3_Figure_3b.eps','epsc')

[r_squared_one_third, p_one_third] = corr(.27 * cf_tab,bdw_tab)
[r_squared_erb, p_erb] = corr(24.7 * (4.37 *  cf_tab/ 1000 + 1),bdw_tab)
[r_squared_fit, p_fit] = corr(fit_.a * log10(cf_tab)+fit_.b,bdw_tab)



%%
cf_tabLOG = log10(cf_tab) ;
bdw_tabLOG = log10(bdw_tab) ;
pow_tabLOG = pow_tab ;

cfTab_vecLOG = log10(cfTab_vec) ;
bdwTab_meanLOG = log10(bdwTab_mean) ;

idxNan = isinf(cfTab_vecLOG) ;
cfTab_vecLOG(idxNan) = [] ;
bdwTab_meanLOG(idxNan) = [] ;

idxNan = isinf(bdwTab_meanLOG) ;
cfTab_vecLOG(idxNan) = [] ;
bdwTab_meanLOG(idxNan) = [] ;

idxNan = cf_tab>1500 ;
bdw_tabLOG(idxNan) = [] ;
cf_tabLOG(idxNan) = [] ;
pow_tabLOG(idxNan) = [] ;

idxNan = isnan(cf_tabLOG) ;
bdw_tabLOG(idxNan) = [] ;
cf_tabLOG(idxNan) = [] ;
pow_tabLOG(idxNan) = [] ;
idxNan = isnan(bdw_tabLOG) ;
bdw_tabLOG(idxNan) = [] ;
cf_tabLOG(idxNan) = [] ;
pow_tabLOG(idxNan) = [] ;

%%
ft = fittype(@(a, x) a + x) ;
fit(1.33*cfTab_vecLOG',bdwTab_meanLOG',ft)

scatter(cf_tabLOG,bdw_tabLOG,.3,'filled');hold on;

onethird_linear = 1.33 * cfTab_vecLOG - 2.39;
plot(cfTab_vecLOG, onethird_linear,'lineWidth',2);


%%
save('./data/results_efficientCoding_all_speech.mat') ;
% %% can be copy and past before figure generation to avoid re-generating all the data
% load('./data/results_efficientCoding_all_speech.mat')
% 
% %%
% cfTab_vec = 2.^(3:.1:13) ;
% 
% powTab_mean = zeros(1,length(cfTab_vec)) ;
% powTab_std  = zeros(1,length(cfTab_vec)) ;
% powTab_nb   = zeros(1,length(cfTab_vec)) ;
% 
% for i = 1:(length(cfTab_vec)-1)
%     ind = find(cf_tab < cfTab_vec(i+1) & cf_tab > cfTab_vec(i)) ;
%     powTab_mean(i) = mean(pow_tab(ind)) ;
%     powTab_std(i)  = std(pow_tab(ind)) ;
%     powTab_nb(i)   = length(pow_tab(ind)) ;
% end
% 
% scatter(cf_tab,pow_tab,1);hold on; 
% 
% set(gca,'xscale','log')
% 
% powTab_mean(isnan(powTab_mean)) = 0 ;
% plot(cfTab_vec,(powTab_mean),'lineWidth',2)
% xlabel('center frequency')
% ylabel('bandwidth')
% legend({'EMD_{raw}','ERB','EMD_{average}'})
% grid on
% hold on;

