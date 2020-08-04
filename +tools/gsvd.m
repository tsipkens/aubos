function [U,V,X,C,S] = gsvd(A,B,flag)
%GSVD   Generalized Singular Value Decomposition.
%   [U,V,X,C,S] = GSVD(A,B) returns unitary matrices U and V,
%   a (usually) square matrix X, and nonnegative diagonal matrices
%   C and S so that
%
%       A = U*C*X'
%       B = V*S*X'
%       C'*C + S'*S = I
%
%   A and B must have the same number of columns, but may have
%   different numbers of rows.  If A is m-by-p and B is n-by-p, then
%   U is m-by-m, V is n-by-n and X is p-by-q where q = min(m+n,p).
%
%   SIGMA = GSVD(A,B) returns the vector of generalized singular
%   values, sqrt(diag(C'*C)./diag(S'*S)).
%
%   The nonzero elements of S are always on its main diagonal.  If
%   m >= p the nonzero elements of C are also on its main diagonal.
%   But if m < p, the nonzero diagonal of C is diag(C,p-m).  This
%   allows the diagonal elements to be ordered so that the generalized
%   singular values are nondecreasing.
%
%   GSVD(A,B,0), with three input arguments and either m or n >= p,
%   produces the "economy-sized" decomposition where the resulting
%   U and V have at most p columns, and C and S have at most p rows.
%   The generalized singular values are diag(C)./diag(S).
%
%   When I = eye(size(A)), the generalized singular values, gsvd(A,I),
%   are equal to the ordinary singular values, svd(A), but they are
%   sorted in the opposite order.  Their reciprocals are gsvd(I,A).
%
%   In this formulation of the GSVD, no assumptions are made about the
%   individual ranks of A or B.  The matrix X has full rank if and only
%   if the matrix [A; B] has full rank.  In fact, svd(X) and cond(X) are
%   are equal to svd([A; B]) and cond([A; B]).  Other formulations, eg.
%   G. Golub and C. Van Loan, "Matrix Computations", require that null(A)
%   and null(B) do not overlap and replace X by inv(X) or inv(X').
%   Note, however, that when null(A) and null(B) do overlap, the nonzero
%   elements of C and S are not uniquely determined.
%
%   Class support for inputs A,B:
%      float: double, single
%
%   See also SVD.

%   Copyright 1984-2016 The MathWorks, Inc.

[m,p]  = size(A);
[n,pb] = size(B);
if pb ~= p
    error(message('MATLAB:gsvd:MatrixColMismatch'))
end

useQA = false;
useQB = false;
if nargin > 2
    if ~(isnumeric(flag) && isscalar(flag) && flag == 0)
        error(message('MATLAB:gsvd:InvalidFlag'));
    end
    % Economy-sized.
    useQA = m > p;
    useQB = n > p;
    if useQA
        [QA,A] = qr(A,0);
        m = p;
    end
    if useQB
        [QB,B] = qr(B,0);
        n = p;
    end
end

[Q,R] = qr([A;B],0);
[U,V,Z,C,S] = csd(Q(1:m,:),Q(m+1:m+n,:));

if nargout < 2
    % Vector of generalized singular values.
    q = min(m+n,p);
    U = [zeros(q-m,1,superiorfloat(A,B)); diagk(C,max(0,q-m))]./[diagk(S,0); zeros(q-n,1)];
else
    % Full composition.
    X = R'*Z;
    if useQA
        U = QA*U;
    end
    if useQB
        V = QB*V;
    end
end


% ------------------------

function [U,V,Z,C,S] = csd(Q1,Q2)
% CSD  Cosine-Sine Decomposition
% [U,V,Z,C,S] = csd(Q1,Q2)
%
% Given Q1 and Q2 such that Q1'*Q1 + Q2'*Q2 = I, the
% C-S Decomposition is a joint factorization of the form
%    Q1 = U*C*Z' and Q2=V*S*Z'
% where U, V, and Z are orthogonal matrices and C and S
% are diagonal matrices (not necessarily square) satisfying
%    C'*C + S'*S = I

[m,p] = size(Q1);
n = size(Q2,1);

if m < n
    [V,U,Z,S,C] = csd(Q2,Q1);
    j = p:-1:1; C = C(:,j); S = S(:,j); Z = Z(:,j);
    m = min(m,p); i = m:-1:1; C(1:m,:) = C(i,:); U(:,1:m) = U(:,i);
    n = min(n,p); i = n:-1:1; S(1:n,:) = S(i,:); V(:,1:n) = V(:,i);
    return
end
% Henceforth, n <= m.

[U,C,Z] = svd(Q1);

q = min(m,p);
i = 1:q;
j = q:-1:1;
C(i,i) = C(j,j);
U(:,i) = U(:,j);
Z(:,i) = Z(:,j);
S = Q2*Z;

if q == 1
    k = 0;
elseif m < p
    k = n;
else
    k = max([0; find(diag(C) <= 1/sqrt(2))]);
end
[V,~] = qr(S(:,1:k));
S = V'*S;
r = min(k,m);
S(:,1:r) = diagf(S(:,1:r));
if m == 1 && p > 1
    S(1,1) = 0;
end

if k < min(n,p)
    r = min(n,p);
    i = k+1:n;
    j = k+1:r;
    [UT,ST,VT] = svd(S(i,j));
    if k > 0
        S(1:k,j) = 0;
    end
    S(i,j) = ST;
    C(:,j) = C(:,j)*VT;
    V(:,i) = V(:,i)*UT;
    Z(:,j) = Z(:,j)*VT;
    i = k+1:q;
    [Q,R] = qr(C(i,j));
    C(i,j) = diagf(R);
    U(:,i) = U(:,i)*Q;
end

if m < p
    % Diagonalize final block of S and permute blocks.
    q = min([nnz(abs(diagk(C,0))>10*m*eps(class(C))), ...
        nnz(abs(diagk(S,0))>10*n*eps(class(C))), ...
        nnz(max(abs(S(:,m+1:p)),[],2)<sqrt(eps(class(C))))]);
    
    % maxq: maximum size of q such that the expression used later on,
    %        i = [q+1:q+p-m, 1:q, q+p-m+1:n],
    % is still a valid permutation.
    maxq = m+n-p;
    q = q + nnz(max(abs(S(:,q+1:maxq)),[],1)>sqrt(eps(class(C))));
    
    i = q+1:n;
    j = m+1:p;
    % At this point, S(i,j) should have orthogonal columns and the
    % elements of S(:,q+1:p) outside of S(i,j) should be negligible.
    [Q,R] = qr(S(i,j));
    S(:,q+1:p) = 0;
    S(i,j) = diagf(R);
    V(:,i) = V(:,i)*Q;
    if n > 1
        i = [q+1:q+p-m, 1:q, q+p-m+1:n];
    else
        i = 1;
    end
    j = [m+1:p 1:m];
    C = C(:,j);
    S = S(i,j);
    Z = Z(:,j);
    V = V(:,i);
end

if n < p
    % Final block of S is negligible.
    S(:,n+1:p) = 0;
end

% Make sure C and S are real and positive.
[U,C] = diagp(U,C,max(0,p-m));
C = real(C);
[V,S] = diagp(V,S,0);
S = real(S);

% ------------------------

function D = diagk(X,k)
% DIAGK  K-th matrix diagonal.
% DIAGK(X,k) is the k-th diagonal of X, even if X is a vector.
if ~isvector(X)
    D = diag(X,k);
    D = D(:);  %Ensure column vector is returned for empty X.
else
    if ~isempty(X) && 0 <= k && 1+k <= size(X,2)
        D = X(1+k);
    elseif ~isempty(X) && k < 0 && 1-k <= size(X,1)
        D = X(1-k);
    else
        D = zeros(0,1,'like',X);
    end
end

% ------------------------

function X = diagf(X)
% DIAGF  Diagonal force.
% X = DIAGF(X) zeros all the elements off the main diagonal of X.
X = triu(tril(X));

% ------------------------

function [Y,X] = diagp(Y,X,k)
% DIAGP  Diagonal positive.
% [Y,X] = diagp(Y,X,k) scales the columns of Y and the rows of X by
% unimodular factors to make the k-th diagonal of X real and positive.
D = diagk(X,k);
j = find(real(D) < 0 | imag(D) ~= 0);
D = diag(conj(D(j))./abs(D(j)));
Y(:,j) = Y(:,j)*D';
X(j,:) = D*X(j,:);
X = X+0; % use "+0" to set possible -0 elements to 0