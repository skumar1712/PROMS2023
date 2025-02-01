% -----------------------------------------------------------------
% In this program, we calculate the singularity spectrum of several
% multifractal signals 
% -----------------------------------------------------------------
% Usage 
% [ ] = XXX()
% Input
%   N:        size of the time series
%   Q:        list of exponents
%   R:
%   wavelet:  analysing wavelet, i.e., 'Gauss', 'DerGauss', 'Sombrero'
%   s0:       minimum scale of analysis (default = 2) 
%   nvoice:   number of voices per octave (default = 10) 
%   noct:     number of octaves discarded by the analysis (default = 2) 
%
% Output
%   h:        Singularity Exponents
%   Dh:       Spectrum of Singularity D(h) calculated using direct method
%             (Lengendre transform)
%   tauq:     Multifractal Exponents tau(q)
%   Z:        Partition function z(q,a)
% -----------------------------------------------------------------

choice = 0 ;
clear
% Load the input file
load('/home/sandekum/Documents/MATLAB/Wavelab850/Books/WaveTour/WTCh06/MATfiles/hdelta_alpha0.9_2N12_qFarey_dt_eq.mat')
X = H_deltat.' ;

N = length(X) ;
Q = linspace(-5, 5, 101) ;
nvoice = 12 ;
s0 = 4 ; c0 = s0; 
noct = 2 ; 
noctave = floor(log2(N)) - noct ;                                          % number of octaves to be considered for analysis 
% wavelet = 'Gauss' ; 
% wavelet = 'DerGauss' ; 
wavelet = 'Sombrero'; 
% wavelet = 'Morlet'; 
X(N+1 : 2*N) = flip(X) ; 
noctave = noctave + 1 ; 
nscale = nvoice * noctave ; 

% Create the file to save the results
% myfile = ['SS_alpha' num2str(1/(1+delta)) 'nvoice' num2str(nvoice) 'N' num2str(log2(N)) '_5Q5.mat'] ;
myfile = ['SS_M' num2str(M) 'nvoice' num2str(nvoice) '2N' num2str(log2(N)) '_5Q5' wavelet '.mat'] ;
% myfile = ['SS_RNDFnvoice' num2str(nvoice) '2N' num2str(log2(N)) '_5Q5.mat'] ;


% The Continuous Wavelet Transform of the signal X
rwt = RWT(X, nvoice, wavelet, noct, s0) ;

% Resize the WT coefficients
[~, nscale] = size(rwt); 
rwt = rwt(:,c0*nvoice + 1 : nscale - nvoice) ;
[n, nscale] = size(rwt) ; 

s = 2.^(log2(s0) + (noctave-c0-1)*(1:nscale)/nscale) ;                     % scales
s = 2.^(1 + (noctave-c0-1)*(1:nscale)/nscale) ;                            % scales

% rwt is arranged so that 1st column correspond to the large scale and last
% column to the small scale 
% noctave-c0-1:  c0 or s0 octaves removed from large scales and 1 from smaller
%                scales 
% Maxima map relating the positions of the WTMM
par = 1e5 ;
maxmap = MM_RWT(rwt, par) ;

% Calculate Thermodynamic parition function 
Z = CalcThermoPartition(rwt, maxmap, Q) ; 
minscale = s(1) ;
maxscale = s(nscale) ;

% Calculate scaling exponent 
tau = FracScalExp(Z, s, minscale, maxscale) ; 
tau = -tau ; 

% Calculate Singularity spectrum 

% hmin = 0; hmax = 1+delta ;                                                 % For Hdelta
% hmin = 0; hmax = 1.5 ;                                                     % For Hdelta
% hmin = 0.4; hmax = 0.8 ;                                                 % For RNDF
hmin = -0.4; hmax = 0.8 ;                                                 % For RNDF
h = hmin : 0.0025 : hmax ;
% h = linspace(hmin, hmax, 101) ; 
Dh = FracSingSpect(tau, Q.', h) ; 
i = find(Dh >= 0);
Dh = Dh(i);    h = h(i);

% Dh_th = h/(1+delta) ;                                                      % For Hdelta
Dh_th = 4*h-2 ;                                                          % For RNDF

Error1 = mean(Dh - Dh_th.') / (hmax - hmin) ;
Error = mean(Dh - Dh_th.') / (hmax - hmin) 
figure; plot(h,Dh,'.')
hold on;
plot(h, Dh_th,'.')

save(myfile)
