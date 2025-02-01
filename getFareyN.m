% ---------------------------------------------------------------------
%                       Calculates the Farey sequence 
% --------------------------------------------------------------------
% Usage 
%   pq = getFareyN(N)
% Input
%   N:      An integer greater than 2
% Ouput
% pq:       A 2D array containing the numerator p and denomenator q 
%           with p, q coprime and q < N.
%           Its length ~ 3 * N^2 / pi^2
% ---------------------------------------------------------------------
% See also
%   gcd, unique, sort
% ---------------------------------------------------------------------
% N = 2^9*3+2; 
N = 117; 
ind = 1 ; 
pq = ones(N^2,2) ; 

if N == 1, error('N should be greater than or equal to 2!'); end
if N == 2, pq = [0 1; 1 1]; end

for qq = 2 : N-1 
    pp = 0 : qq   ;
%     pq = zeros(N+1,2); 
    for j = 1 : length(pp)
        if gcd(pp(j),qq) == 1
            p = pp(j)  ;
            q = qq  ;
        else
            p = pp(j) / gcd(pp(j),qq) ; 
            q = qq / gcd(pp(j),qq) ; 
        end
        pq(ind,:) =[p q] ;
        ind = ind+1 ;
    end
     
end
pq_aux = unique(pq,'rows');
[s,ii] = sort(pq_aux(:,1) ./ pq_aux(:,2)); 
pq = pq_aux(ii,:) ; 

myfile = ['gcdpqN' num2str(N) '.mat' ] ;
save(myfile, 'N','pq') ;
 
% end