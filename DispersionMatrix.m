function [DMat] = DispersionMatrix(n, dispRate)
%DISPERSIONMATRIX Builds Dispersion matrix for given n and dispRate.
%   Uses tensor products to create a sparse n^2 x n^2 matrix D such that
%   D*C (where each column of C represents the contents of an nxn grid)
%   gives C after one tick of dispersion.


isq2 = 1/sqrt(2);
%sets up matrices to handle dispersion within the same row and with
%adjacent row
%sameRowMat = gallery('tridiag', n, dispRate, 1-(4+4*isq2)*dispRate, dispRate);
%offRowMat = gallery('tridiag', n, isq2*dispRate, dispRate, isq2*dispRate);
sameRowMat = toeplitz([1-(4+4*isq2)*dispRate, dispRate, zeros(1,n-3), dispRate]);
offRowMat = toeplitz([dispRate, isq2*dispRate, zeros(1, n-3), isq2*dispRate]);
%indicators for where those matrices go
sameRowDiag = spdiags(ones(n,1), 0, n, n);
offRowDiag = spdiags(ones(n-1,1), -1, n, n);
%assembling the dispersion matrix with tensor products
sameDisp = kron(sameRowDiag, sameRowMat);
downDisp = kron(offRowDiag, offRowMat);
DMat = downDisp + sameDisp + downDisp';
%%figure(3); spy(DMat)
%corrections for the corners and edges to prevent dispersion out of bounds
cornersMat = sparse(n,n);
edgesMat = speye(n);
cornersMat(1,1) = 1;
cornersMat(n,n) = 1;
edgesMat(1,1) = 0;
edgesMat(n,n) = 0;
cornersBig = kron(cornersMat, cornersMat);
edgesBig = kron(edgesMat, cornersMat) + kron(cornersMat, edgesMat);
%%figure(4); spy(cornersBig);
%%figure(5); spy(edgesBig);
%DMat = DMat + edgesBig*(1+2*isq2)*dispRate + cornersBig*(2+3*isq2)*dispRate;
DMat = DMat + kron(cornersMat, edgesMat)*(1+2*isq2)*dispRate;

end

