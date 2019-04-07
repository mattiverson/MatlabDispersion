%Problem Parameters
R1 = [1 1; 0 2];
R2 = [0 2; 1 1];
k1 = 3.5e-5;
k2 = 1.0e-8;

%Computation Parameters
n = 20;
ticksToRun = 1000;
dispRate = 0.05;
isq2 = 1/sqrt(2);

%Animation Parameters
framerate = 5;
period = 20;

%Initial State
range = linspace(0,1,n);
[x, y] = meshgrid(range, range);


%Dispersion Matrix (Don't touch)

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

%Tracking which reactions depend on which chemicals (Don't touch)
r1Dep = zeros(1,3);
r2Dep = zeros(1,3);
r1Dep((R1(1,1) ~= 0) + 2*(R1(1,2) ~= 0)) = 1;
r2Dep((R2(1,1) ~= 0) + 2*(R2(1,2) ~= 0)) = 1;

C = zeros(n^2,2);
C(1,1) = n^2 / 2;
C(end,2) = n^2 / 2;
shading interp

t = 0;
startTime = clock;
while t<ticksToRun
    %dispersion
    C = DMat*C;
    
    %reaction 1
    if r1Dep(1)
        r1Rate = C(:,1)*k1;
        r1In = kron(r1Rate, R1(1,:));
        r1Out = kron(r1Rate, R1(2,:));
    elseif r1Dep(2)
        r1Rate = C(:,2)*k1;
        r1In = kron(r1Rate, R1(1,:));
        r1Out = kron(r1Rate, R1(2,:));
    else
        r1Rate = C(:,1).*C(:,2)*k1;
        r1In = kron(r1Rate, R1(1,:));
        r1Out = kron(r1Rate, R1(2,:));
    end
    
    %reaction 2
    if r2Dep(1)
        r2Rate = C(:,1)*k2;
        r2In = kron(r2Rate, R2(1,:));
        r2Out = kron(r2Rate, R2(2,:));
    elseif r2Dep(2)
        r2Rate = C(:,2)*k2;
        r2In = kron(r2Rate, R2(1,:));
        r2Out = kron(r2Rate, R2(2,:));
    else
        r2Rate = C(:,1).*C(:,2)*k2;
        r2In = kron(r2Rate, R2(1,:));
        r2Out = kron(r2Rate, R2(2,:));
    end
    
    %add back reaction results
    C = C - r1In - r2In + r1Out + r2Out;
    t = t + 1;
    %draw results
    if mod(t, period) == 0
        fprintf("Tick %d done \n", t);
        figure(1)
        surf(x,y, reshape(log(C(:,1) + 1e-12),n,n))
        shading interp;
        view(2)
        %caxis([0 1]);
        colorbar;
        title("Log of Concentration of Chemical 1")
        figure(2)
        surf(x,y, reshape(log(C(:,2) + 1e-12),n,n))
        shading interp;
        view(2)
        %caxis([0 1]);
        title("Log of Concentration of Chemical 2")
        colorbar;
        figure(3)
        surf(x,y, reshape(log(r1Rate + 1e-12),n,n))
        shading interp;
        view(2)
        title("Log of Rate of Reaction 1")
        colorbar;
        figure(4)
        surf(x,y, reshape(log(r2Rate + 1e-12),n,n))
        shading interp;
        view(2)
        title("Log of Rate of Reaction 2")
        colorbar;
        endTime = clock;
        duration = endTime(6)-startTime(6);
        fprintf("%.1f", duration*1000);
        pause(max(1/framerate - duration, 0));
        %pause;
        startTime = clock;
    end
end
disp("Done")