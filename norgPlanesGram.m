function G = norgPlanesGram(A, eps)
%   norgPlanesGram Compute canonical Gramian between nowhere orthogonal
%   Euclidean planes
%   
%   INPUT:
%   A: Gramian of orthobases of n nowhere orthogonal planes in R^d
%   eps: [optional] Tolerance for testing of input
%
%   OUTPUT:
%   G: Gramian of orthobases of the same n nowhere orthogonal planes in R^d 
%   which is in a form which is an injective invariant for nowhere 
%   orthogonal tuples in the Grassmannian of planes in R^d modulo O(d)
%
%   CITATION:
%   Algorithm 1 in
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


% Testing Input
n=round(size(A,1)/2);

if nargin == 1
    eps = 1e-10;  % To do: optimize the default tol
end

if ~isreal(A) || (size(A,1) ~= size(A,2)) || (size(A,1) ~= 2*n) || ...
        (norm(A-A') > eps) || (trace(kron(eye(n),[0 0; 1 0])*A) > eps) ||...
        (abs(max(diag(A))-mean(diag(A))) > eps)
    error('Input must be the Gramian of orthobases of n planes in R^d')
end

R = [1 0; 0 -1];

H = zeros(size(A));

distsvFlag=0;

for ii=1:n
    for jj=ii:n
        [W, Sig, V]=svd(A((2*(ii-1)+1):2*ii,(2*(jj-1)+1):2*jj));
        if Sig(1,1)*Sig(2,2) < eps
            error('Input must be the Gramian of orthobases of nowhere orthogonal planes in R^d')
        elseif abs(Sig(1,1)-Sig(2,2)) < eps
            H((2*(ii-1)+1):2*ii,(2*(jj-1)+1):2*jj)=W*V';
        else 
            distsvFlag=1;
            
            H((2*(ii-1)+1):2*ii,(2*(jj-1)+1):2*jj)=W*V';
            
            break
        end
    end
    if distsvFlag
        break
    end
end

H=H+(H-diag(diag(H)))';

if distsvFlag %  there exists (kk,ll) s.t. A_{kk,ll} has distinct sing. val.
    kk=ii; 
    ll=jj;
    
    D=zeros(size(A));
    Dtil=zeros(size(A));
    
    [Wk, ~, ~]=svd(A((2*(ii-1)+1):2*ii,(2*(jj-1)+1):2*jj));
    Wktil=Wk*R;
        
    for jj=1:(kk-1)
        D((2*(jj-1)+1):2*jj,(2*(jj-1)+1):2*jj)=...
            H((2*(kk-1)+1):2*kk,(2*(jj-1)+1):2*jj)'*Wk;
        Dtil((2*(jj-1)+1):2*jj,(2*(jj-1)+1):2*jj)=...
            H((2*(kk-1)+1):2*kk,(2*(jj-1)+1):2*jj)'*Wktil;
    end
    
    D((2*(kk-1)+1):2*kk,(2*(kk-1)+1):2*kk)=Wk;
    Dtil((2*(kk-1)+1):2*kk,(2*(kk-1)+1):2*kk)=Wktil;
    
    for jj=(kk+1):ll
        D((2*(jj-1)+1):2*jj,(2*(jj-1)+1):2*jj)=...
            H((2*(kk-1)+1):2*kk,(2*(jj-1)+1):2*jj)'*Wk;
        Dtil((2*(jj-1)+1):2*jj,(2*(jj-1)+1):2*jj)=...
            H((2*(kk-1)+1):2*kk,(2*(jj-1)+1):2*jj)'*Wktil;
    end
    
    for jj=(ll+1):n
        [W, ~, V]=svd(A((2*(kk-1)+1):2*kk,(2*(jj-1)+1):2*jj));
        Hkj=W*V';
        D((2*(jj-1)+1):2*jj,(2*(jj-1)+1):2*jj)=Hkj'*Wk;
        Dtil((2*(jj-1)+1):2*jj,(2*(jj-1)+1):2*jj)=Hkj'*Wktil;
    end
    
    G1=D'*A*D;
    G2=Dtil'*A*Dtil;
else
    negdetFlag=0;
    
    D=zeros(size(A));
    
    S=kron(eye(n),R);
    
    for jj=1:n
        D((2*(jj-1)+1):2*jj,(2*(jj-1)+1):2*jj)=...
            H(1:2,(2*(jj-1)+1):2*jj);
    end
    
    DHD=D*H*D';
    
    for ii=1:n
        for jj=ii:n
            if det(DHD((2*(ii-1)+1):2*ii,(2*(jj-1)+1):2*jj))<0 % there exists (kk,ll) s.t det((DHD*)_{kk,ll})=-1
                negdetFlag=1;
                
                kk=ii;
                ll=jj;
                
                Q=inv(sqrtm(DHD((2*(kk-1)+1):2*kk,(2*(ll-1)+1):2*ll)*R));
                E=kron(eye(n),Q)*D;
                
                G1=E*A*E';
                G2=S*E*A*E'*S;
                
                break
            end
        end
        if negdetFlag
            break
        end
    end
    
    if ~negdetFlag
        G1=D*A*D';
        G2=S*D*A*D'*S;
    end
end

% returning (row-dominant) lexicographic min
diffG=G1'-G2';
if diffG(find(abs(diffG)>eps,1))<0
    G=G1;
else
    G=G2;
end

end

