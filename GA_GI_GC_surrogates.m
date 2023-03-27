function out=GA_GI_GC_surrogates(S,t,d,p,Fs,band_all,numsurr,perc_GC_GI,perc_GA,q,nfft)

%%% INPUT
% S --- time series NxM (N = length; M = no. of series)
% t --- target series
% d --- driver series
% p --- model order subject-specific
% Fs --- sampling freq subject-specific
% band_all --- frequency bands of interest - each band is in a row
%               if one does not want to compute freq. measures, specify band_all = 0
% numsurr --- no. of surrogates
% perc --- percentiles [low median high]
% q --- no. of lags for correlation functions
% nfft --- no. of points (positive frequencies)

%%% OUTPUT
% out structure with prctiles of F_XY and A_Y (time and frequency measures)

%% check input
if nargin < 11, nfft = 1000; end
if nargin < 10, q = 20; end
if nargin < 9, perc_GA = [2.5 50 97.5]; end
if nargin < 8, perc_GC_GI = [5 50 95]; end
if nargin < 7, numsurr = 100; end
if nargin < 6, error('Please specify frequency band to compute surrogates'); end
if nargin < 5, error('Please specify sampling frequency'); end
if nargin < 4, error('Please specify model order'); end
if nargin < 3, error('Not enough parameters'); end

%% init
b=size(band_all,1);

A_Y_surr=zeros(1,numsurr);
F_XY_surr=A_Y_surr;
F_Y_surr=A_Y_surr;

A_Y_band_surr=zeros(b,numsurr);
F_XY_band_surr=A_Y_band_surr;
F_Y_band_surr=A_Y_band_surr;

f_XY_surr=zeros(nfft,numsurr);
a_Y_surr=f_XY_surr;
f_Y_surr=f_XY_surr;

A_Y_band_s=zeros(b,length(perc_GA));
F_XY_band_s=A_Y_band_s;
F_Y_band_s=A_Y_band_s;

%% VAR model identification
[~,M]=size(S);

% restricted: driver from target & driver
out_dr=LinReg(S,d,[1 2],(1:p));
eAm_dr=out_dr.eA';
%eSu_dr=out_dr.es2u;
U_dr=out_dr.eu;

% restricted for GA: target from driver
% --> we are neglecting the internal dynamics of target
out_GA_targ=LinReg(S,t,d,(1:p)); % maintain the same model order 
eAm_targ_GA=out_GA_targ.eA;
%eSu_targ_GA=out_GA_targ.es2u;
U_targ_GA=out_GA_targ.eu;

%%% put coefficients eAm_dr & eAm_targ_GA in one coefficient matrix eAm_GA
eAm_GA=zeros(M,M*p); 
eAm_GA(d,:)=eAm_dr;
if t==1
    for ord_tmp=1:p
        eAm_GA(t,2*ord_tmp)=eAm_targ_GA(ord_tmp);
    end
    %eSu_GA=[eSu_targ_GA 0; 0 eSu_dr];
     U_GA=[U_targ_GA U_dr];
else
    for ord_tmp=1:p
        eAm_GA(t,2*ord_tmp-1)=eAm_targ_GA(ord_tmp);
    end
    %eSu_GA=[eSu_dr 0; 0 eSu_targ_GA];
    U_GA=[U_dr U_targ_GA];
end

% restricted for GC: target from target
% ---> we are neglecting the causal coupling from driver to target
out_GC_targ=LinReg(S,t,t,(1:p));
eAm_targ_GC=out_GC_targ.eA;
%eSu_targ_GC=out_GC_targ.es2u;
U_targ_GC=out_GC_targ.eu;

%%% put coefficients eAm_dr & eAm_targ_GC in one coefficient matrix eAm_GC:
eAm_GC=zeros(M,M*p);
eAm_GC(d,:)=eAm_dr;
if t==1
    for ord_tmp=1:p
        eAm_GC(t,2*ord_tmp-1)=eAm_targ_GC(ord_tmp);
    end
    %eSu_GC=[eSu_targ_GC 0; 0 eSu_dr];
    U_GC=[U_targ_GC U_dr];
else
    for ord_tmp=1:p
        eAm_GC(t,2*ord_tmp)=eAm_targ_GC(ord_tmp);
    end
    %eSu_GC=[eSu_dr 0; 0 eSu_targ_GC];
    U_GC=[U_dr U_targ_GC];
end

%% cycle on all surrogates
for k = 1:numsurr
    % resampling of the AR residuals 
    % P = randperm(N) returns a vector containing a random permutation of the
    % integers 1:N.  For example, randperm(6) might be [2 4 5 6 1 3].
    l=length(U_GA);
    idx = randperm(l);
    U_surr_GA=U_GA(idx,:); % same resampling for the two series
    U_surr_GC=U_GC(idx,:); 

    %%% GA surrogates
    S_surr_GA=MVARfilter(eAm_GA,U_surr_GA'); % same estimated AR coefficients, resampled residuals
    out_surro_GA_full=LinReg(S_surr_GA',[1 2],[1 2],(1:p)); % VAR model identification (full)
    eAm_GAsurr=out_surro_GA_full.eA';
    eSu_GAsurr=out_surro_GA_full.es2u;
    ret_GAsurr = GC_GI_GA_computation(eAm_GAsurr,eSu_GAsurr,t,d,q,nfft,Fs); % compute measures for the surrogate    
    A_Y_surr(k)=ret_GAsurr.A_Y; %%% GA time domain
    a_Y_surr(:,k)=ret_GAsurr.a_Y_all;  %%% GA frequency domain
    
    %%% GC surrogates
    S_surr_GC=MVARfilter(eAm_GC,U_surr_GC'); % same estimated AR coefficients, resampled residuals
    out_surro_GC_full=LinReg(S_surr_GC',[1 2],[1 2],(1:p)); % VAR model identification (full)
    eAm_GCsurr=out_surro_GC_full.eA';
    eSu_GCsurr=out_surro_GC_full.es2u;
    ret_GCsurr = GC_GI_GA_computation(eAm_GCsurr,eSu_GCsurr,t,d,q,nfft,Fs); % compute measures for the surrogate
    F_XY_surr(k)=ret_GCsurr.F_XY; %%% GC time domain
    f_XY_surr(:,k)=ret_GCsurr.f_XY; %%% GC frequency domain
    %%% GI surrogates
    F_Y_surr(k)=ret_GCsurr.F_Y; %%% GI time domain
    f_Y_surr(:,k)=ret_GCsurr.f_Y; %%% GI frequency domain

    if band_all ~= 0
        %%% band values:
        for i=1:b
            band=band_all(i,:);
            A_Y_band_surr(i,k)=sum(a_Y_surr(band(1):band(2),k))/nfft;
            F_XY_band_surr(i,k)=sum(f_XY_surr(band(1):band(2),k))/nfft;
            F_Y_band_surr(i,k)=sum(f_Y_surr(band(1):band(2),k))/nfft;
        end
    end
    
end

%%% confidence intervals
F_XY_s = prctile(F_XY_surr,perc_GC_GI);
F_Y_s = prctile(F_Y_surr,perc_GC_GI);
A_Y_s = prctile(A_Y_surr,perc_GA);
a_Y_s=prctile(a_Y_surr,perc_GA,2); % freq. distribution of GA
f_XY_s = prctile(f_XY_surr,perc_GC_GI,2); % freq. distribution of GC
f_Y_s=prctile(f_Y_surr,perc_GC_GI,2); % freq. distribution of GI

for i=1:b
    tmp_A_Y=A_Y_band_surr(i,:);
    tmp_F_XY=F_XY_band_surr(i,:);
    tmp_F_Y=F_Y_band_surr(i,:);
    A_Y_band_s(i,:) = prctile(tmp_A_Y,perc_GA);
    F_XY_band_s(i,:) = prctile(tmp_F_XY,perc_GC_GI);
    F_Y_band_s(i,:) = prctile(tmp_F_Y,perc_GC_GI);
end

% return values in output
%%%%%% GC
out.F_XY_s=F_XY_s;
out.F_XY_band_s=F_XY_band_s;
out.f_XY_s=f_XY_s;
%%%%%% GA
out.A_Y_s=A_Y_s;
out.A_Y_band_s=A_Y_band_s;
out.a_Y_s=a_Y_s;
%%%%%% GI
out.F_Y_s=F_Y_s;
out.F_Y_band_s=F_Y_band_s;
out.f_Y_s=f_Y_s;

end
