%% APPLICATION TO CEREBROVASCULAR INTERACTIONS MAP-CBFV
clear; close all; clc;

%% Data analysis parameters
numsurr=100; % number of surrogates
perc_GC_GI=[5 50 95]; % percentiles for GC and GI surrogates
perc_GA=[2.5 50 97.5]; % percentiles for GA surrogates

nfft=1000; % number of points on frequency axis 
pfilter=0.94; % filter parameter for detrending
q=20; % number of lags for computing autocorrelation functions
p=7; % model order (fixed for the representative example)
rangeLF1=[0.02 0.07]; % VLF band
rangeLF2=[0.07 0.2]; % LF band

idata=[3 4]; % 3) MAP; 4) CBFV 
t=2; % target index
d=1; % driver index

% plot parameters
col1=[192 0 0]./255;
col2=[0 80 150]./255;
DimensioneFont=17;
axislinewidth=0.5;

%% load data
load('data.mat');
So=data(:,idata); % extract the series
M=size(So,2); % this is bivariate analysis

% sampling frequency as the inverse of the heart period HP
HP=data(:,1); 
fs=1/mean(HP);
% select spectral band of interest
band1=round((nfft*2/fs)*rangeLF1);
band2=round((nfft*2/fs)*rangeLF2);

% pre-processing: AR filtering and removal of mean value
Sf=nan*ones(size(So,1),M);
S=Sf;
for m=1:M
    Sf(:,m)=GICA_ARfilter(So,m,pfilter); % AR highpass filtered series
    S(:,m)=Sf(:,m)-mean(Sf(:,m)); % zero-mean series
end
            
%%  IDENTIFICATION OF THE FULL LINEAR MODEL
jv=[1 2]; % index of predicted series 
iv=[1 2]; % indexes of predictors 
iv_lags=(1:p); % lags of predictors
outARX=GICA_LinReg(S,jv,iv,iv_lags);
Am=outARX.eA; Su=outARX.es2u; % estimated ARX parameters

%% get GC, GI and GA measures in time and frequency domains
ret = GICA_computation(Am',Su,t,d,q,nfft,fs);
f=ret.f; % frequency axis

F_XY=ret.F_XY; % Granger Causality GC
f_XY=ret.f_XY; % frequency GC
F_Y=ret.F_Y; % Granger Isolation (GI)
f_Y=ret.f_Y; % frequency GI
A_Y=ret.A_Y; % Granger autonomy (GA)
a_Y=ret.a_Y_all; % frequency GA 

% average values in bands
F_XY_band1=sum(f_XY(band1(1):band1(2)))/nfft;
F_XY_band2=sum(f_XY(band2(1):band2(2)))/nfft;
F_Y_band1=sum(f_Y(band1(1):band1(2)))/nfft;
F_Y_band2=sum(f_Y(band2(1):band2(2)))/nfft;
A_Y_band1=sum(a_Y(band1(1):band1(2)))/nfft;
A_Y_band2=sum(a_Y(band2(1):band2(2)))/nfft;

% Decomposition of spectrum of the full model
P=ret.P; % spectral matrix
P_t = abs(squeeze(P(t,t,:))); % target spectrum 
P_d = abs(squeeze(P(d,d,:))); % driver spectrum

%% surrogate data analysis
band_f=[band1; band2];
out=GICA_surrogates(S,t,d,p,fs,band_f,numsurr,perc_GC_GI,perc_GA,q,nfft);

F_XY_s=out.F_XY_s; % time GC
A_Y_s=out.A_Y_s; % time GA
F_Y_s=out.F_Y_s; % time GI
a_Y_s=out.a_Y_s; % this is the vector of percentiles (for each freq. point)
f_XY_s=out.f_XY_s; % this is the vector of percentiles (for each freq. point)
f_Y_s=out.f_Y_s; % this is the vector of percentiles (for each freq. point)

F_XY_band_s=out.F_XY_band_s; % freq GC 
A_Y_band_s=out.A_Y_band_s; % freq GA 
F_Y_band_s=out.F_Y_band_s; % freq GI 

% exctract information in band:
A_Y_band1_s=A_Y_band_s(1,:); A_Y_band2_s=A_Y_band_s(2,:);
F_XY_band1_s=F_XY_band_s(1,:); F_XY_band2_s=F_XY_band_s(2,:);
F_Y_band1_s=F_Y_band_s(1,:); F_Y_band2_s=F_Y_band_s(2,:);

% count number of significant values
%%%% granger autonomy
if A_Y > A_Y_s(3); A_Y_count_surr=1;
else; A_Y_count_surr=0; end
if A_Y_band1 > A_Y_band1_s(3) || ...
        abs(A_Y_band1) > abs(A_Y_band1_s(1))
    A_Y_band1_count_surr=1;
else; A_Y_band1_count_surr=0;
end
if A_Y_band2 > A_Y_band2_s(3) || ...
        abs(A_Y_band2) > abs(A_Y_band2_s(1))
    A_Y_band2_count_surr=1;
else; A_Y_band2_count_surr=0;
end
%%% granger causality
if F_XY > F_XY_s(3); F_XY_count_surr=1;
else; F_XY_count_surr=0; end
if F_XY_band1 > F_XY_band1_s(3); F_XY_band1_count_surr=1;
else; F_XY_band1_count_surr=0; end
if F_XY_band2 > F_XY_band2_s(3); F_XY_band2_count_surr=1;
else; F_XY_band2_count_surr=0; end
%%% granger isolation
if F_Y < F_Y_s(1); F_Y_count_surr=1;
else; F_Y_count_surr=0; end
if F_Y_band1 < F_Y_band1_s(1); A_Y_band1_count_surr=1;
else; F_Y_band1_count_surr=0; end
if F_Y_band2 < F_Y_band2_s(1); F_Y_band2_count_surr=1;
else; F_Y_band2_count_surr=0; end

%% plot

h0=figure('numbertitle','off',...
    'WindowState','maximized','Color','w');
figure(h0); 

spl{1}=subplot(2,5,[1 2]); % MAP
plot(Sf(:,d),'color',col2,'LineWidth',1.5)
axis([0 length(Sf(:,d)) 0.98*min(Sf(:,d)) 1.02*max(Sf(:,d))]);
ylabel('MAP {\it [mmHg]}'); xticks([]);
legend('X_n'); legend box off
set(gca,'FontSize',DimensioneFont,'FontName','Times');
ax=gca;
ax.LineWidth=axislinewidth;
set(spl{1}, 'Position', [0.06,0.53,0.25,0.38]);
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k','linewidth',axislinewidth,'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k','linewidth',axislinewidth,'HandleVisibility','off')
if isholdonque == 0
hold off
end
text(0,ax.YLim(2)+((ax.YLim(2)-ax.YLim(1))/10),...
    'A','FontWeight','bold','FontSize',DimensioneFont+2)

spl{2}=subplot(2,5,[6 7]); % CBFV
plot(Sf(:,t),'color',col1,'LineWidth',1.5)
axis([0 length(Sf(:,t)) 0.98*min(Sf(:,t)) 1.02*max(Sf(:,t))]);
ylabel('MCBFV {\it [cm/s]}');  xlabel('{\it n} [beats]');
legend('Y_n'); legend box off
box on;
set(gca,'FontSize',DimensioneFont,'FontName','Times');
ax=gca;
ax.LineWidth=axislinewidth;
set(spl{2}, 'Position', [0.06,0.11,0.25,0.38]);
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k','linewidth',axislinewidth,'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k','linewidth',axislinewidth,'HandleVisibility','off')
if isholdonque == 0
hold off
end

spl{3}=subplot(2,5,3); % MAP spectrum
plot(f,P_d,'color',col2,'LineWidth',2); 
xticks(0:0.1:fs/2); xticks([]);
ymax=1.01*max(P_d);
axis([0 0.5 0 ymax]);
ylabel('PSD_{MAP}{\it [mmHg^2]}');
legend('P_X');
legend('boxoff')
set(gca,'FontSize',DimensioneFont,'FontName','Times');
ax=gca;
ax.LineWidth=axislinewidth;
set(spl{3}, 'Position', [0.37,0.53,0.15,0.38]);
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k','linewidth',axislinewidth,'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k','linewidth',axislinewidth,'HandleVisibility','off')
if isholdonque == 0
hold off
end
text(0,ax.YLim(2)+((ax.YLim(2)-ax.YLim(1))/10),...
    'B','FontWeight','bold','FontSize',DimensioneFont+2)

spl{4}=subplot(2,5,8); % CBFV spectrum
plot(f,P_t,'color',col1,'LineWidth',2); hold on;
xticks([0 0.1 0.2 0.3 0.4]); xlabel('{\it f}[Hz]');
ymax=1.01*max(P_t);
axis([0 0.5 0 ymax]);
ylabel('PSD_{MCBFV}{\it [(cm/s)^2]}');
legend('P_Y');
legend('boxoff')
set(gca,'FontSize',DimensioneFont,'FontName','Times');
ax=gca;
ax.LineWidth=axislinewidth;
set(spl{4}, 'Position', [0.37,0.11,0.15,0.38]);
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k','linewidth',axislinewidth,'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k','linewidth',axislinewidth,'HandleVisibility','off')
if isholdonque == 0
hold off
end

spl{5}=subplot(2,5,4); % GC
area(f,f_XY_s(:,3),'FaceColor',[0.9 0.9 0.9], ...
    'HandleVisibility','off'); hold on;
area(f,f_XY_s(:,1),'FaceColor',[1 1 1], ...
    'HandleVisibility','off'); hold on;
plot(f,f_XY_s(:,1),'LineWidth',1,'Color',[0.5 0.5 0.5],...
    'HandleVisibility','off'); hold on; % low perc
plot(f,f_XY_s(:,3),'LineWidth',1,'Color',[0.5 0.5 0.5],...
    'HandleVisibility','off'); hold on; % high perc
plot(f,f_XY_s(:,2),'LineWidth',1.5,'Color',[0.2 0.2 0.2],...
    'HandleVisibility','off'); hold on;  % 50° perc
plot(f,f_XY,'Color',[0.4660 0.6740 0.1880], ...
    'LineWidth',2); hold on; 
xticks([]); 
xlim([0 0.5]);
yticks([0 1 2 3])
legend('{\it f}_{X\rightarrowY}');
legend('boxoff')
set(gca,'FontSize',DimensioneFont,'FontName','Times');
ax=gca;
ax.LineWidth=axislinewidth;
set(spl{5}, 'Position', [0.565,0.53,0.18,0.38]);
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k','linewidth',axislinewidth,'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k','linewidth',axislinewidth,'HandleVisibility','off')
if isholdonque == 0
hold off
end
text(0,ax.YLim(2)+((ax.YLim(2)-ax.YLim(1))/10),...
    'C','FontWeight','bold','FontSize',DimensioneFont+2)

spl{6}=subplot(2,5,9); % GI
area(f,f_Y_s(:,3),'FaceColor',[0.9 0.9 0.9], ...
    'HandleVisibility','off'); hold on;
area(f,f_Y_s(:,1),'FaceColor',[1 1 1], ...
    'HandleVisibility','off'); hold on;
plot(f,f_Y_s(:,1),'LineWidth',1,'Color',[0.5 0.5 0.5], ...
    'HandleVisibility','off'); hold on; % low perc
plot(f,f_Y_s(:,3),'LineWidth',1,'Color',[0.5 0.5 0.5], ...
    'HandleVisibility','off'); hold on; % high perc
plot(f,f_Y_s(:,2),'LineWidth',1.5,'Color',[0.2 0.2 0.2], ...
    'HandleVisibility','off'); hold on; % 50° perc
plot(f,f_Y,'-','Color',[0.3010 0.7450 0.9330], ...
    'LineWidth',2); hold on; 
xlabel('{\it f}[Hz]');
xticks(0:0.1:fs/2);
xlim([0 0.5]);
ylim([0 11]);
yticks([0 3 6 9])
legend('{\it f}_{Y}')
legend('boxoff')
set(gca,'FontSize',DimensioneFont,'FontName','Times');
ax=gca;
ax.LineWidth=axislinewidth;
set(spl{6}, 'Position', [0.565,0.11,0.18,0.38]);
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k','linewidth',axislinewidth,'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k','linewidth',axislinewidth,'HandleVisibility','off')
if isholdonque == 0
hold off
end

spl{7}=subplot(2,5,[5 10]); % GA
area(f,a_Y_s(:,3),'FaceColor',[0.9 0.9 0.9], ...
    'HandleVisibility','off'); hold on;
area(f,a_Y_s(:,1),'FaceColor',[0.9 0.9 0.9], ...
    'HandleVisibility','off'); hold on;
plot(f,a_Y_s(:,1),'LineWidth',1,'Color',[0.5 0.5 0.5], ...
    'HandleVisibility','off'); hold on; % low perc
plot(f,a_Y_s(:,3),'LineWidth',1,'Color',[0.5 0.5 0.5], ...
    'HandleVisibility','off'); hold on; % high perc
plot(f,a_Y_s(:,2),'LineWidth',1.5,'Color',[0.2 0.2 0.2], ...
    'HandleVisibility','off'); hold on; % 50° perc
plot(f,a_Y,'Color',[0.9290 0.6940 0.1250],'LineWidth',2);  hold on;
xlabel('{\it f}[Hz]');
xticks(0:0.1:fs/2);
xlim([0 0.5]);
ylim([-0.9 3]);
yticks([0 1 2 3])
legend('{\it a}_{Y}')
legend('boxoff')
set(gca,'FontSize',DimensioneFont,'FontName','Times');
ax=gca;
ax.LineWidth=axislinewidth;
set(spl{7}, 'Position', [0.79,0.11,0.18,0.8]);
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k','linewidth',axislinewidth,'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k','linewidth',axislinewidth,'HandleVisibility','off')
if isholdonque == 0
hold off
end
text(0,ax.YLim(2)+((ax.YLim(2)-ax.YLim(1))/25),...
    'D','FontWeight','bold','FontSize',DimensioneFont+2)

