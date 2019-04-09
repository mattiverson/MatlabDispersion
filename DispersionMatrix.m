function [DMat] = DispersionMatrix(n, dispRate)
%DISPERSIONMATRIX Summary of this function goes here
%   Detailed explanation goes here

%sets up matrices to handle dispersion within the same row and with
%adjacent row
sameRowMat = gallery('tridiag', n, dispRate, 1-(4+4*isq2)*dispRate, dispRate);
offRowMat = gallery('tridiag', n, isq2*dispRate, dispRate, isq2*dispRate);
%indicators for where those matrices go
sameRowDiag = spdiags(ones(n,1), 0, n, n);
offRowDiag = spdiags(ones(n-1,1), -1, n, n);
%assembling the dispersion matrix with tensor products
sameDisp = kron(sameRowDiag, sameRowMat);
downDisp = kron(offRowDiag, offRowMat);
DMat = downDisp + sameDisp + downDisp';
%%figure(3); spy(DMat)
%corrections for the corners and edges to prevent dispersion out of bounds
cornersMat = zeros(n,n);
edgesMat = diag(ones(n,1));
cornersMat(1,1) = 1;
cornersMat(n,n) = 1;
edgesMat(1,1) = 0;
edgesMat(n,n) = 0;
cornersBig = kron(cornersMat, cornersMat);
edgesBig = kron(edgesMat, cornersMat) + kron(cornersMat, edgesMat);
%%figure(4); spy(cornersBig);
%%figure(5); spy(edgesBig);
DMat = DMat + edgesBig*(1+2*isq2)*dispRate + cornersBig*(2+3*isq2)*dispRate;


end

