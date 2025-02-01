% ----------------------------------------------------------------
%                   Compute $H_\delta(t)$ as in KPVV2021
% ----------------------------------------------------------------
% Usage 
% 
% Input 
%   J:              determines the length of the signal N = 2^J
%   alpha:          alpha = 1/(1+delta)
%   gcdpqNn.mat:    Contains the Farey sequence p/q, q<n
% 
% Output
%   hdelta:         1xN dimensional array: h_\delta (Dirac comb with 
%                   delta's located at the Farey sequence)
%   Hdelta:         1xN dimensional array: H_\delta  
% ----------------------------------------------------------------
% Description
%   KPVV2021: Section 4
% ----------------------------------------------------------------
clear; 
tic
% close all; 

for k = 3 : 2 : 9  
    %        load the Farey sequences 
    load('gcdpqN165.mat')
    pq2 = pq ;
    load('gcdpqN164.mat')
    pq1 = pq ;
    
    J = 9; N=2^(J+4) ;                                                     % Length of the signal 
    %       Get the rationals in pq2 but not in pq1 along with their positions 
    %       in the parent set pq2
    [c,ii] = setdiff(pq2(:,1)./pq2(:,2), pq1(:,1)./pq1(:,2)); 
    %       Compute the number of remaining rationals n
    n = 2^(J+4)-size(pq1,1)+1 ;
    a = [0.51, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95 0.99];         % The alpha values,  with 1/2 < alpha < 1
    alpha = a(k) ;

    delta = 1/alpha - 1;
    n = 1 ; 
    b_ndelta = newgamma((n+2*delta)/2) ...
           / ((2*pi^(n/2+2*delta)) *abs(newgamma(-delta)))   ;             % newgamma calculates the Gamma function 


    %       Choose the additional n rationals randomly from the set c 
    iin = randperm(numel(ii),n) ; 
    pq1_n = [pq ; pq2(ii(iin),1) pq2(ii(iin),2)] ; 
    [s, id] = sort(pq1_n(:,1) ./ pq1_n(:,2)) ; 
    p = pq1_n(id,1) ;
    q = pq1_n(id,2) ;

    %       Calculate the coefficients of the Dirac comb
    coeff = zeros(length(p),1) ;
    norm2_zeta = 1/sqrt(2) ; 
    for j = 1 : length(p)
        if mod(q(j),2) == 1
            coeff(j) = 1 / q(j)^(2*(1+delta)) ;
        elseif mod(q(j)/2,2) == 1 
            coeff(j) = -2 * (2^(1+2*delta)-1) / q(j)^(2*(1+delta)) ;
        else
            coeff(j) = 2^(2*(1+delta)) / q(j)^(2*(1+delta)) ;
        end
    end 
    %       Calculate h_pdelta
    h_pdelta = -2 *b_ndelta * zeta(2+2*delta) * coeff / norm2_zeta; 

    t1 = linspace(0,1,N+1) ;  
    H_deltat = zeros(1,N);
    %       Calculate H_delta 
    for j = 1 : N
        ss1 = s(s<=t1(j)) ;                                                    % all the rational p/q less than t(j)
        [ii1, id1] = ismember(ss1, s) ;                                        % look for their location, i.e., id in the parent array 's'
        H_deltat(j) = sum(h_pdelta(id1)) ;                                     % and sum those h_pdelta's         
    end

    Hdeltat1 = cumsum(h_pdelta) ;                                              % Calculating by summing all the deltas
    
    myfile = ['hdelta_alpha' num2str(alpha) '_2N' num2str(log2(N)) '_qFarey_dt_eq.mat' ] ;
    save(myfile, 'h_pdelta','H_deltat', 'Hdeltat1', 'p', 'q', 's','delta') ;

end
toc   

% ----------------------------------------------------------------
% Date:             February 15, 2021
% Updated:          Dec 10, 2021
% Written by:       Sandeep Kumar
% ----------------------------------------------------------------
