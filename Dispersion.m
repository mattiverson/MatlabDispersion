%Problem Parameters
R(:,:,1) = [1 1; 0 2];
R(:,:,2) = [0 2; 1 1];
k = [3.5e-2 1e-3];

%Computation Parameters
n = 40;
ticksToRun = 150;
dispRate = 0.05;

%Animation Parameters
period = ticksToRun - 1;

%Initial State
range = linspace(0,1,n);
[x, y] = meshgrid(range, range);
C1Initial = 2*sin(pi*x) .* sin(pi*y);
C2Initial = sin(pi*x);
C = zeros(n^2, 2);
C(:,1) = reshape(C1Initial, n^2, 1);
C(:,2) = reshape(C2Initial, n^2, 1);

%Dispersion Matrix 
DMat = DispersionMatrix(n,dispRate);

%Tracking which reactions depend on which chemicals
rDep = R(1,:,:) ~= 0;

t = 0;
while t<ticksToRun
    
    %dispersion
    C = DMat*C;
    
    %reactions
    rIn = zeros(size(C));
    rOut = zeros(size(C));
    for(idx = 1:length(R(1,1,:)))
        rRate = prod(C(:,rDep(:,:,idx)), 2) * k(idx);
        rIn = rIn + kron(rRate, R(1,:,idx));
        rOut = rOut + kron(rRate, R(2,:,idx));
    end

    C = C - rIn + rOut;
    
    %draw results
    %performace note: the graphics code takes roughly 1000x longer than
    %the computations in a tick
    if mod(t, period) == 0
        figure(1)
        %colormap gray
        surf(x,y, reshape(C(:,1),n,n))
        shading interp;
        view(2)
        title("Log of Concentration of Chemical 1")
        colorbar;
        caxis([0 1]);
        
        figure(2)
        %colormap gray
        surf(x,y, reshape(C(:,2),n,n))
        shading interp;
        view(2)
        title("Log of Concentration of Chemical 2")
        colorbar;
        caxis([0 1]);
    end
    t = t + 1;
end
