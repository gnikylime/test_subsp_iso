function areIso = traceWords(A,B,W,V,eps)
%   traceWords 
%   For n-tuples of (d x d) matrices (A_1,...,A_n) and (B_1,...,B_n), 
%   determine if there exists a unitary U s.t. U A_i U' = B_i for all i.
%       
%   INPUT:
%   A: d x (dn) matrix representing an n-tuple of (d x d) matrices
%   B: d x (dn) matrix representing an n-tuple of (d x d) matrices to be
%       compared to A.
%   W: Cell array of words in 1, 2, ..., n (represented as vectors) s.t. 
%       the words evaluated on the A_1, A_2, ..., A_n is a basis for a
%       desired algebra A((A_i)_i).  E.g., the output of canBasMatAlg.
%   V: Cell array of words in 1, 2, ..., n (represented as vectors) s.t. 
%       the words evaluated on the B_1, B_2, ..., B_n is a basis for a
%       desired algebra A((B_i)_i).  Must be generated in the same order as W.
%       E.g., the output of canBasMatAlg.
%   eps: [optional] Tolerance for testing of equality of trace
%
%   OUTPUT:
%   areIso: Bool indicating whether the subspace tuples are
%       
%   CITATION:
%   Lemma 11 in
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
%   TO DO:
%   Improve testing of input W for elements being vectors.
%   Improve testing of A, B to ensure closure under adjoints.

% Testing Input
d = size(A,1);
n = round(size(A,2)/d);

if (size(A,2) ~= d*n) || (size(B,1) ~= d) || (size(B,2) ~= d*n)
    error('First two inputs must d x (dn) matrices.')
end

if ~iscell(W)
    error('The third input must be a cell array of numerical vectors')
end

if nargin == 4
    eps = 1e-10;  % To do: optimize the default tol
end

% Lemma 11
m=length(W);

E = [eye(d) zeros(d,d*(m-1))]; % word evals on A matrices
F = [eye(d) zeros(d,d*(m-1))]; % word evals on B matrices

for jj=2:m
    word = W{jj};
    Ejj = eye(d);
    Fjj = eye(d);
    for ll = 1:length(word)
        Ejj = Ejj * A(:,(d*(ll-1)+1):d*ll);
        Fjj = Fjj * B(:,(d*(ll-1)+1):d*ll);
    end
    E(:,(d*(jj-1)+1):d*jj) = Ejj;
    F(:,(d*(jj-1)+1):d*jj) = Fjj;
end

if ~isequal(W,V)
    disp('Condition (i) failed.')
    areIso = 0;
    return;
end

for ii = 1:m
    for jj = 1:m
        Eij = E(:,(d*(ii-1)+1):d*ii)'*E(:,(d*(jj-1)+1):d*jj);
        Fij = F(:,(d*(ii-1)+1):d*ii)'*F(:,(d*(jj-1)+1):d*jj);
        if (abs(trace(Eij-Fij))> eps)
            disp('Condition (ii) failed.')
            areIso = 0;
            return
        end
        for kk = 1:m
            if (abs(trace(Eij*E(:,(d*(kk-1)+1):d*kk)-...
                    Fij*F(:,(d*(kk-1)+1):d*kk))) > eps)
                disp('Condition (iii) failed.')
                areIso = 0;
                return;
            end
        end
    end
    for ll = 1:n
        if (abs(trace(E(:,(d*(ii-1)+1):d*ii)'*A(:,(d*(ll-1)+1):d*ll)...
                -F(:,(d*(ii-1)+1):d*ii)'*B(:,(d*(ll-1)+1):d*ll))) > eps)
            disp('Condition (iv) failed.')
            areIso = 0;
            return;
        end
    end
end

areIso = 1;

end