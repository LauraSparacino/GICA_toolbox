function ret = GC_GI_GA_computation(Am,Su,t,d,q,nfft,Fs)

%   --- PART 1 --- PART 2 ---
%   Compute Information Dynamics analytically for a stationary 2AR(p) process:
%   S_n=A_1*S_{n-1}+A_2*S_{t-n}+...+A_p*S_{n-p}+E_n (**)
%   where S = vector containing two series, d & t (X & Y, one is the driver and the other is the target)
%   __ Instead of using the data, uses the 2AR parameters (either estimated from data, or theoretical)
%   __ Works under the linear Gaussian assumption

%   --- PART 3 --- PART 4 --- PART 5 ---
%   Compute spectral Granger autonomy (GA), Granger causality (GC) and
%   Granger Isolation (GI)
%   starting from the spectral representation of the above 2AR model (full) 
%   and the spectral representation of the following 2AR model containing
%   a restricted regression of Y on X (X model):
%   X_n = A_xx(z)*X_{n-k} + A_xy(z)*Y_{n-k} + U_x|xy (*)
%   Y_n = B_yx(z)*X_{n-k} + ............... + U_y|x
%   In this restricted model, we assume S=[X Y], i.e. X-driver is
%   indicated by index 1 and Y-target by index 2

%   INPUT: 
%   Am  -   generalized connectivity matrix A=(A_1 A_2 ... A_p)
%   Su  -   covariance matrix for E_n
%   t   -   target series
%   d   -   driver series
%   q   -   number of lags (at which you want to compute correlations - q>p)
%   nfft - number of points for calculation of the spectral functions (nfft)
%   Fs - sampling frequency
%   OUTPUT:
%   ret structure with all possible conditional entropies,
%   frequency-specific Granger autonomy and Granger causality

M = size(Am,1); % number of elements in system
assert(M==2); % this is bivariate analysis, error if M~=2
p=floor(size(Am,2)/M); % number of lags in 2AR model (**)

if nargin<5, Fs = 1; end
if nargin<4, nfft = 1000; end     
if all(size(nfft)==1) % if nfft is scalar
    f = (0:nfft-1)*(Fs/(2*nfft)); % frequency axis
else % if nfft is a vector, we assume that it is the vector of the frequencies
    f = nfft; nfft = length(nfft);
end
% zeta is the variable on the complex plane
z = 1i*2*pi/Fs;

%% init
A = [eye(M) -Am]; % build (I-Am)
G=zeros(M,M,nfft); % Transfer Matrix of (*)
P_r=zeros(M,M,nfft); % Spectral Matrix of (*)
H=zeros(M,M,nfft); % Transfer Matrix of (**)
P_f=zeros(M,M,nfft); % Spectral Matrix of (**)
caus_r=zeros(nfft,1); caus_f=caus_r; % causal spectra
aut_f=caus_r; aut_r=caus_r; % autonomous spectra
ga_freq_all=caus_r; % spectral GA - sum of two terms
ga_freq_variable=caus_r; % spectral GA - term variable in frequency
gc_freq=caus_r; % spectral GC
gi_freq=caus_r; % spectral GI

%% PART 1: COVARIANCE MATRICES
% prepare covariance matrices, (:,:,1) is lag 0, (:,:,q+1) is lag q 
R=NaN*ones(M,M,q+1); 

% Obtain F (actually A^p) and Delta (actually Sigma^p)
Im=eye(M*p);
F=[Am;Im(1:end-M,:)];% this is A^p
Delta=zeros(p*M,p*M); % this is actually Sigma^p, but use Delta for clarity in code
Delta(1:M,1:M)=Su(:,:);

% Obtain R_o^p=BigSigma solving the
% Lyapunov equation: BigSigma = F * BigSigma * F^T + Delta
BigSigma=dlyap(F,Delta);

% extract R(0),...,R(p-1) from R_o^p
% R(0) contains the variances on the main diagonal
for i=1:p
    R(:,:,i)=BigSigma(1:M,M*(i-1)+1:M*i);
end

% Yule-Walker solution for lags >= p (until lag q)
for k=p+1:q+1
    Rk=R(:,:,k-1:-1:k-p);
    Rm=[];
    for ki=1:p
        Rm=[Rm; Rk(:,:,ki)];
    end
    R(:,:,k)=Am*Rm;
end

L=q; % maximum lag of Xpast, Ypast

%%% create covariance matrices related to X,Y
Ry=NaN*ones(1,q+1);
Rx=NaN*ones(1,q+1);
Rxy=NaN*ones(1,q+1);
Ryx=NaN*ones(1,q+1);
for k=1:q+1
    Ry(k)=R(t,t,k); % target
    Rx(k)=R(d,d,k); % driver
    Rxy(k)=R(d,t,k); % X given Y
    Ryx(k)=R(t,d,k); % Y given X
end

%% PART 2: CONDITIONAL ENTROPIES given Y as target
%%% Entropy of Y
Sy=Ry(1);
Hy=0.5*log(Sy)+0.5*log(2*pi*exp(1));

%%% Conditional Entropy of Y given Ypast - AR model
% T = toeplitz(c,r) returns a nonsymmetric Toeplitz matrix with c as its first column
% and r as its first row. 
SYpast=toeplitz(Ry(1:L),Ry(1:L));
SyYpast=Ry(2:L+1);
Sy_Ypast = Sy - SyYpast/SYpast * SyYpast';
Hy_Ypast=0.5*log(Sy_Ypast)+0.5*log(2*pi*exp(1));

%%% Conditional Entropy of Y given Xpast - X model
SXpast=toeplitz(Rx(1:L),Rx(1:L));
SyXpast=Ryx(2:L+1);
Sy_Xpast = Sy - SyXpast/SXpast * SyXpast';
Hy_Xpast=0.5*log(Sy_Xpast)+0.5*log(2*pi*exp(1));

%%% Conditional Entropy of Y given Ypast and Xpast - ARX model
SYX=toeplitz(Rxy(1:L),Ryx(1:L));
SXY=SYX';
SYpastXpast=[[SYpast SYX];[SXY SXpast]];
SyYpastXpast=[SyYpast SyXpast];
Sy_YpastXpast = Sy - SyYpastXpast/SYpastXpast * SyYpastXpast';
Hy_YpastXpast=0.5*log(Sy_YpastXpast)+0.5*log(2*pi*exp(1));

%% PART 2: CONDITIONAL ENTROPIES given X as target
%%% Entropy of X
Sx=Rx(1);
Hx=0.5*log(Sx)+0.5*log(2*pi*exp(1));

%%% Conditional Entropy of X given Xpast - AR model
SXpast=toeplitz(Rx(1:L),Rx(1:L));
SxXpast=Rx(2:L+1);
Sx_Xpast = Sx - SxXpast/SXpast * SxXpast'; 
Hx_Xpast=0.5*log(Sx_Xpast)+0.5*log(2*pi*exp(1));

%%% Conditional Entropy of X given Ypast - X model
SYpast=toeplitz(Ry(1:L),Ry(1:L));
SxYpast=Rxy(2:L+1);
Sx_Ypast = Sx - SxYpast/SYpast * SxYpast';
Hx_Ypast=0.5*log(Sx_Ypast)+0.5*log(2*pi*exp(1));

%%% Conditional Entropy of X given Ypast and Xpast - ARX model
SXY=toeplitz(Ryx(1:L),Rxy(1:L));
SYX=SXY';
SXpastYpast=[[SXpast SXY];[SYX SYpast]];
SxXpastYpast=[SxXpast SxYpast];
Sx_XpastYpast = Sx - SxXpastYpast/SXpastYpast * SxXpastYpast';
Hx_XpastYpast=0.5*log(Sx_XpastYpast)+0.5*log(2*pi*exp(1));

%% PART 3: RESTRICTED MODEL of Y (X model) & FULL MODEL of X (ARX model) --- (*)
% extract coefficients and partial variances related to the model (*):
%   X_n = A_xx(z)*X_{n-k} + A_xy(z)*Y_{n-k} + U_x|xy
%   Y_n = B_yx(z)*X_{n-k} + ............... + U_y|x
    
%%%%% UPDATE indices of target & driver to refer to (*)
t_restricted=2; d_restricted=1;

B = SyXpast / SXpast ;
B=B(1:p); % coefficients X --> Y (X model)
% Sy_Xpast -- partial variance of X model (X-->Y) 
% Sx_XpastYpast -- partial variance of ARX model (X,Y-->X) 
S_tilde_r=[Sx_XpastYpast 0; 0 Sy_Xpast];

for n=1:nfft 
    %%% Coefficient matrix in the frequency domain
    As = zeros(M,M); % matrix As(z)=I-sum(A(k))
    Bs = 0; % matrix Bs(z)=sum(B(k))
    for k = 1:p+1
        As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n));
        % indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix B
        % (A(1) is in the second block, and so on)
    end
    for k = 1:p
        Bs = Bs + B(k)*exp(-z*(k)*f(n));
    end

    AB=[As(d,[d t]); -Bs 1]; % transfer matrix for restricted model
    G(:,:,n)  = inv(AB);
    P_r(:,:,n)  = G(:,:,n)*S_tilde_r*G(:,:,n)'; 

    % "causal" part of the spectrum P_r:
    caus_r(n)=Sx_XpastYpast*abs(G(t_restricted,d_restricted,n))^2;
    % "autonomous" part of the spectrum P_r:
    aut_r(n)=Sy_Xpast*abs(G(t_restricted,t_restricted,n))^2;

end

%% PART 4: FULL MODEL (2AR) (**)
% extract partial variances
% Sy_YpastXpast -- partial variance of 2AR model (X,Y-->Y)
% Sx_XpastYpast -- partial variance of 2AR model (X,Y-->X)
if t==1
    S_tilde_f=[Sy_YpastXpast 0; 0 Sx_XpastYpast];
else
    S_tilde_f=[Sx_XpastYpast 0; 0 Sy_YpastXpast];
end

for n=1:nfft 
    %%% Coefficient matrix in the frequency domain
    As = zeros(M,M); % matrix As(z)=I-sum(A(k))
    for k = 1:p+1
        As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n));
        % indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix B
        % (A(1) is in the second block, and so on)
    end

    H(:,:,n)  = inv(As);
    P_f(:,:,n)  = H(:,:,n)*S_tilde_f*H(:,:,n)'; 

    % "causal" part of the spectrum P_f:
    caus_f(n)=Sx_XpastYpast*abs(H(t,d,n))^2;
    % "autonomous" part of the spectrum P_f:
    aut_f(n)=Sy_YpastXpast*abs(H(t,t,n))^2;
end

%% PART 5: GRANGER AUTONOMY & CAUSALITY
for n=1:nfft

    % SPECTRAL GC --- Geweke formulation (without 1/2)
    gc_freq(n)=log(abs(P_f(t,t,n))./aut_f(n)); 
    gi_freq(n)=log(abs(P_f(t,t,n))./caus_f(n));

    % SPECTRAL GA
    num=Sy_Xpast*abs(H(t,t,n))^2;
    den=Sy_YpastXpast*abs(G(t_restricted,t_restricted,n))^2;
    ga_freq_all(n)=log(num./den); 
    ga_freq_variable(n)=log((abs(H(t,t,n))^2)./(abs(G(t_restricted,t_restricted,n))^2));

end

% time domain measures
gc_time=sum(gc_freq)/nfft; %GC
gi_time=sum(gi_freq)/nfft; %GI
ga_time=sum(ga_freq_all)/(nfft); %GA

%% OUTPUTs
ret.f=f;
ret.P=P_f;

%%%% target Y
ret.Hy=Hy; ret.Sy=Sy;
ret.Hy_y=Hy_Ypast; ret.Sy_y=Sy_Ypast;
ret.Hy_x=Hy_Xpast; ret.Sy_x=Sy_Xpast;
ret.Hy_yx=Hy_YpastXpast; ret.Sy_yx=Sy_YpastXpast; 

%%%% target X
ret.Hx=Hx; ret.Sx=Sx;
ret.Hx_x=Hx_Xpast; ret.Sx_x=Sx_Xpast;
ret.Hx_y=Hx_Ypast; ret.Sx_y=Sx_Ypast;
ret.Hx_xy=Hx_XpastYpast; ret.Sx_xy=Sx_XpastYpast; 

%%%% partial spectra of (*) and (**) models
ret.caus_r=caus_r;
ret.aut_r=aut_r;
ret.caus_f=caus_f;
ret.aut_f=aut_f;

%%%% GRANGER CAUSALITY & AUTONOMY
% freq
ret.a_Y_all=ga_freq_all;
ret.a_Y_variable=ga_freq_variable;
ret.f_XY=gc_freq;
ret.f_Y=gi_freq;
% time
ret.A_Y=ga_time;
ret.F_XY=gc_time;
ret.F_Y=gi_time;

%%%% COEFFICIENTS OF THE RESTRICTED MODEL Y_X and TRANSFER FUNCTIONS
ret.G=G; 
ret.H=H; 
ret.B=B;

end