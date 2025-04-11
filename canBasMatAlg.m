function W = canBasMatAlg(A)
%   canBasMatAlg 
%   Canonical basis for matrix algebra from finite generating set
%   
%   INPUT:
%   A: k x (kn) matrix where the, A((A_i)_i) is the smallest algebra with
%       unity that contains the n (k x k) submatrices of A
%
%   OUTPUT:
%   W: Cell array of words in 1, 2, ..., n (represented as vectors) s.t. 
%       the words evaluated on the A_1, A_2, ..., A_n is a basis for A((A_i)_i)
%       
%   CITATION:
%   Algorithm 2 in
%   Testing isomorphism between tuples of subspaces
%   Emily J King, Dustin G Mixon, Shayne Waldron
%   2021
%   arXiv:2105.03448
%
%   CONTACT:
%   Emily J King, emily.king@colostate.edu
%
%   LAST UPDATED:
%   November 10, 2021
%
%   TO DO LIST:
%   Decide if a more numerically stable test than rank is desired in
%   practice.

% Testing Input
k = size(A,1);
n = round(size(A,2)/k);

if (size(A,2) ~= k*n) 
    error('Input must a k x (kn) matrix')
end

% Implementation of Algorithm 2
W = {[]}; % the word of length zero
wordEval = reshape(eye(k),[],1); % evaluation of word of length zero

m_old = 0;
m_new = 1;

while m_new > m_old
    m_old = m_new;
    for ii = 1:n
        for jj = 1:m_old
            tmp = reshape(A(:,(k*(ii-1)+1):k*ii)*reshape(wordEval(:,jj),k,k),[],1);
            if rank(wordEval) ~= rank([wordEval tmp])
                W{m_new+1} = [ii W{jj}];
                wordEval(:,m_new+1) = tmp;
                m_new=m_new+1;  
            end
        end
    end
end

end
