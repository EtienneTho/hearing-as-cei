% Etienne Thoret 2023 (c)
%
% Generate figure 1b: Separation of tonal signal mixtures by Cascaded
% Envelope Interpolation
%
% If you use this script please cite the following paper
% Thoret, E., Ystad, S., Kronland-Martinet, R. (2023) Hearing as adaptive cascaded envelope interpolation
% Communications Biology
%

close all ;
clearvars;
addpath(genpath('./lib/'));

%%
N = 16000;
t = 1:N;

p = N/5;

% chirp
f0 = .01 ;
phi = rand ;
x2 = amgauss(N,N/2,N/4).*chirp(t,f0,16000,f0+.18,'linear',phi)' ;

% fmod
f0 = f0+.4 ;
x1 = amgauss(N,N/2,N/2).*fmsin(N,f0,f0+.05,p,N/2,f0+.05);

x = real(x1+x2);
x = x;

nbIMFs = 6 ;
nbSiftingIterations = 1 ;
imf = emdc_fix([],x, nbSiftingIterations, nbIMFs) ;            


thresh = -15 ;

figure;
subplot(331)
plot(real(x1),'color','k')
ylabel('Amplitude')
xlabel('Time')
title('x_1')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])
axis('square')
axis([0 N -2 2])


subplot(332)
plot(real(x2),'color','k')
ylabel('Amplitude')
xlabel('Time')
title('x_2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])
axis('square')
axis([0 N -2 2])


subplot(333)
plot(real(x),'color','k')
ylabel('Amplitude')
xlabel('Time')
title('x_1 + x_2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])
axis('square')
axis([0 N -2 2])

subplot(334)
plot(imf(1,:),'color','k');hold on;
%plot(real(x1)'-imf(1,:),'color','r');hold off;
ylabel('Amplitude')
xlabel('Time')
title('Mode 1')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])
axis('square')
axis([0 N -2 2])

subplot(335)
plot(imf(2,:),'color','k');hold on;

ylabel('Amplitude')
xlabel('Time')
title('Mode 2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])
axis('square')
axis([0 N -2 2])

subplot(336)
plot(sum(imf(1,:)+imf(2,:),1),'color','k');hold on;
ylabel('Amplitude')
xlabel('Time')
title('Mode 1 + Mode 2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])
axis('square')
axis([0 N -2 2])

subplot(337)
spectrogram(imf(1,:),256,16,'yaxis','MinThreshold',thresh)
% s_imf1 = spectrogram(imf(1,:),256,16,'yaxis','MinThreshold',thresh);
% imagesc(flipud(abs(s_imf1)));
ylabel('Frequency')
xlabel('Time')
colorbar('off')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Mode 1')
axis('square')

subplot(338)
spectrogram(imf(2,:),256,16,'yaxis','MinThreshold',thresh)
ylabel('Frequency')
xlabel('Time')
colorbar('off')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Mode 2')
axis('square')

subplot(339)
s_imf1 = spectrogram(imf(1,:),256,16,'yaxis','MinThreshold',thresh);
s_imf2 = spectrogram(imf(2,:),256,16,'yaxis','MinThreshold',thresh);
imagesc(flipud(abs(s_imf1+s_imf2)));

ylabel('Frequency')
xlabel('Time')
colorbar('off')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Mode 1 + Mode 2')
axis('square')

set(gcf,'position',[300,300,700,500])

%%
saveas(gca,'F1_Figure_1b.eps','epsc')
saveas(gca,'F1_Figure_1b.fig')
