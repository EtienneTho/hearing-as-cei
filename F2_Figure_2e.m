% Etienne Thoret 2023 (c)
%
% Generate figure 2e: Combination tone, theoretical vs. simulated
% spectrum
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%
% Figure_2e.csv - col 1: Theoretical cubic CT, col 2: CEI cubic CT 


% this script runs the analysis of distorsion products by EMD
close all ;
clearvars;
addpath(genpath('./lib/'));

fs =  16000 ;
duration = .2 ;% # of data samples
Nt = floor(fs*duration) ;
t = linspace(0,duration,Nt) ;

Nf = 5 ;
Nratio = 30 ;

f1 = linspace(1000,4000,Nf) ;
% f1 = 1000
f1 = [500 750 1000 1250 1500 1750 2000] ;
ratio = linspace(1.3, 1.7, Nratio) ;

Nratio = length(ratio) ;
Nf = length(f1) ;

% compute the frequency of the distorsion product from the EMD
% decomposition
Nrepet = 1 ;
f_2f1_f2_emd_tab_tab   = zeros(Nrepet,Nf,Nratio) ;
amp_2f1_f2_emd_tab_tab = zeros(Nrepet,Nf,Nratio) ; 
tab_sum_tab_tab        = zeros(Nrepet,Nf,Nratio) ; 

for iRepet = 1:Nrepet
    iRepet
    f_2f1_f2_tab       = zeros(Nf,Nratio) ;
    ratio_tab          = zeros(Nf,Nratio) ;
    f_2f1_f2_emd_tab   = zeros(Nf,Nratio) ;
    amp_2f1_f2_emd_tab = zeros(Nf,Nratio) ;    
    tab_sum            = zeros(Nf,Nratio) ;    
    for i_f1 = 1:Nf
        for i_ratio = 1:Nratio
            i_f1/Nf ;
            i_ratio/Nratio ;
            f1(i_f1) ;
            f1(i_f1) * ratio(i_ratio) ;
            f_2f1_f2_tab(i_f1,i_ratio) = 2*f1(i_f1) - f1(i_f1) * ratio(i_ratio) ;
            ratio_tab(i_f1,i_ratio) = ratio(i_ratio) ;
            phi1 = rand*2*pi ;
            phi2 = rand*2*pi ;
            s1 = cos(2*pi*f1(i_f1)*t+phi1) ;
            s2 = cos(2*pi*f1(i_f1) * ratio(i_ratio)*t+phi2) ;
            x = (s1 + s2) ;
            [f_2f1_f2_emd_loop,amp_2f1_f2_emd_loop, locations_loop, peaks_loop, f__loop] = distortion_products(x, fs, f1(i_f1)) ;
            f_2f1_f2_emd_tab(i_f1,i_ratio) = f_2f1_f2_emd_loop ;
            amp_2f1_f2_emd_tab(i_f1,i_ratio) = amp_2f1_f2_emd_loop ;
            
            tab_sum(i_f1,i_ratio) = sum(abs (floor(f__loop) - floor(f_2f1_f2_tab(i_f1,i_ratio)))<=1) ;
        end

    end
    f_2f1_f2_emd_tab_tab(iRepet,:,:)   = f_2f1_f2_emd_tab ;
    amp_2f1_f2_emd_tab_tab(iRepet,:,:) = amp_2f1_f2_emd_loop ;
    tab_sum_tab_tab(iRepet,:,:) = tab_sum ;
end

f_2f1_f2_emd_tab   = squeeze(mean(f_2f1_f2_emd_tab_tab,1)) ;
amp_2f1_f2_emd_tab = squeeze(mean(amp_2f1_f2_emd_tab_tab,1)) ;
tab_sum_tab_tab    = squeeze(mean(tab_sum_tab_tab,1)) ;

%%
figure
plot((0:1500),(0:1500))
hold on
scatter(f_2f1_f2_tab(:),f_2f1_f2_emd_tab(:))

addpath('./lib') ;
opt.XLabel = 'Theoretical cubic CT (Hz)'; % xlabel
opt.YLabel = 'CEI cubic CT (Hz)'; %ylabel
opt.BoxDim = [5, 2]; %[width, height]

axis([0 1500 0 1500]);

% apply
setPlotProp(opt);

%%
saveas(gca,'./out/eps/F2_Figure_2e.eps','epsc')
saveas(gca,'./out/fig/F2_Figure_2e.fig','fig')

% Figure_2e.csv - col 1: Theoretical cubic CT, col 2: CEI cubic CT 
matrixtowrite = [f_2f1_f2_tab(:), f_2f1_f2_emd_tab(:)];
writematrix(matrixtowrite, './out/csv/Figure_2e.csv') ;

%%

function [f_2f1_f2_emd, amp_2f1_f2_emd, locations, peaks, f_out] = distortion_products(x, fs, fmin)
    nbIMFs = 6 ;
    nbSiftingIterations = 1 ;
    nbSampleEnvelope = 1 ;
    imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;    
    [imf_emd] = emdc([],x);
    [szRow,~] = size(imf) ;
    p_spec = [] ;
    
    for iImfs = 1:szRow
         h_ = hilbert(imf(iImfs,:)) ;
        [p__, f__] = pspectrum(imf(iImfs,:), fs) ;
        p_spec = [p_spec p__] ;
    end
    
    emdemd_pspect = (mean(pspectrum(imf_emd'),2));

    idxFmin = find(f__<fmin) ;
    idxFmin = max(idxFmin) ;

    emd_pspectrum = (mean(p_spec,2)) ;
    [peaks,locations] = findpeaks(emd_pspectrum(1:idxFmin), 'SortStr', 'descend') ;
    [sorted_peaks, sorted_indexes] = sort(peaks) ;
    
    f_2f1_f2_emd = f__(locations(sorted_indexes(end-2))) ;
    amp_2f1_f2_emd = max(peaks(sorted_indexes(end)),peaks(sorted_indexes(end-1))) - peaks(sorted_indexes(end-2)) ;
    f_2f1_f2_emd = f__(locations(3)) ;
    amp_2f1_f2_emd = peaks(1)-peaks(3) ;
    f_2f1_f2_emd = f__(locations(1)) ;
    amp_2f1_f2_emd = 0 ;
    f_out = f__(locations) ;

end




