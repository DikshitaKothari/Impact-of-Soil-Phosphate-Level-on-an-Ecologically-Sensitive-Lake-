clc  
close all 

plotFlag = 1;

% Spatial (x) parameters analytical
L      = 80000;                        % Size of the domain 
ystart = 0;                         % Start of computational domain (m)
yend   = ystart + L;                % End of computational domain (m)
nya    = 1000;                        % Spatial grid points 
dya     = (yend - ystart)/(nya - 1);  % Spatial step size, delta r (m)
ya      = ystart: dya : yend;         % Vector of grid points

% Spatial y parameters numerical
ny =  [9, 17, 33, 65, 129, 257];
for i = 1:length(ny)
    dy(i) = (yend - ystart)/(ny(i) - 1);  % Spatial step size, delta r (m) 
end
y1 = ystart: dy(1) : yend;
y2 = ystart: dy(2) : yend;
y3 = ystart: dy(3) : yend;
y4 = ystart: dy(4) : yend;
y5 = ystart: dy(5) : yend;
y6 = ystart: dy(6) : yend;

% Phyiscal parameters
phi0     = 50;                          % Temperature at Y=0
phi_1    = 0;                           % Temperature at Y=L
D        = 2.37*10^2;                   % Diffusion equation constant                 
d        = sqrt(D);                     % Square root of Diffusion constant

% Temporal (t) parameters 
timestart = 0;                              % Starting time (0 seconds)
timeend   = 86400*28;                       % Ending Time 28 days
dt        = 206;                            % Temporal step size, delta t (s)
nt        = (timeend - timestart)/dt + 1;   % Time grid points
time      = timestart:dt:timeend;           % Vector of times

% Initialise the solution arrays - Tn (T^n) and Tnp1 (T^(n+1)) and exact
% solution TExact
Cn1 = zeros(1,ny(1));
Cn2 = zeros(1,ny(2));
Cn3 = zeros(1,ny(3));
Cn4 = zeros(1,ny(4));
Cn5 = zeros(1,ny(5));
Cn6 = zeros(1,ny(6));

Cnp11 = zeros(1,ny(1));
Cnp12 = zeros(1,ny(2));
Cnp13 = zeros(1,ny(3));
Cnp14 = zeros(1,ny(4));
Cnp15 = zeros(1,ny(5));
Cnp16 = zeros(1,ny(6));

CExact = zeros(1,nya);    % Exact fourier series solution for temperature
%%
% Specify the initial conditions (already zeros)

% Enforce the boundary conditions - These will never change in the simulation
Cn1(1)       = 50;
Cn2(1)       = 50;
Cn3(1)       = 50;
Cn4(1)       = 50;
Cn5(1)       = 50;
Cn6(1)       = 50;

Cnp11(1)     = 50;
Cnp12(1)     = 50;
Cnp13(1)     = 50;
Cnp14(1)     = 50;
Cnp15(1)     = 50;
Cnp16(1)     = 50;


%Fourier parameters - pre-compute B and lambda vectors
nfs     = 100;                  % Number of fourier terms               
B       = zeros(1,nfs);         % Initialise B vector   
lambda  = zeros(1,nfs);         % Initialise lambda vector
p       = zeros(1,nfs);         % Initialise p vector
for n = 1:nfs;           
    p(n)        = (n-0.5)*pi/L;
    B(n)        = -(2*phi0/(L*p(n)));                   
    lambda(n)   = d*p(n);
end


%% Main solution loop

% Loop through time
for i=2:nt                   % loop starts from 2 (not 1) to ensure that the analytical and numerical solutions are compared at the same time point in the mean error calc below
    
    t=time(i);               % current time for outputted solution
    
    CExact = zeros(1,nya);

    %Analytical solution (same as Week 3)
    for j = 1:nya                  % For each grid point        
        CExact(j) = phi0;         % Add the steady state solution
        
        % Loop through the Fourier series
        for n = 1:nfs
            CExact(j)  = CExact(j) + B(n)*sin(p(n)*ya(j))*exp(-lambda(n)^2*t);  % Calculate series sum for T at r(i) and t
        end
    end
    
    % 1 % Numerical solution - note that we only need to solve for j=2 to nr-1
    % since the first and last points are fixed by the boundary condition
    sigma   =  D*dt/(dy(1)^2);
    for j   =   2:ny(1)-1
        Cnp11(j) =   sigma*Cn1(j+1)+(1-2*sigma)*Cn1(j)+sigma*Cn1(j-1);
    end
        Cnp11(ny(1)) =   Cn1(ny(1))+sigma*(Cn1(ny(1)-1)-Cn1(ny(1)));
    
        % 2 % Numerical solution - note that we only need to solve for j=2 to nr-1
    % since the first and last points are fixed by the boundary condition
    sigma   =  D*dt/(dy(2)^2);
    for j   =   2:ny(2)-1
        Cnp12(j) =   sigma*Cn2(j+1)+(1-2*sigma)*Cn2(j)+sigma*Cn2(j-1);
    end
        Cnp12(ny(2)) =   Cn2(ny(2))+sigma*(Cn2(ny(2)-1)-Cn2(ny(2)));
    
        % 3 % Numerical solution - note that we only need to solve for j=2 to nr-1
    % since the first and last points are fixed by the boundary condition
    sigma   =  D*dt/(dy(3)^2);
    for j   =   2:ny(3)-1
        Cnp13(j) =   sigma*Cn3(j+1)+(1-2*sigma)*Cn3(j)+sigma*Cn3(j-1);
    end
        Cnp13(ny(3)) =   Cn3(ny(3))+sigma*(Cn3(ny(3)-1)-Cn3(ny(3)));
    
        % 1 % Numerical solution - note that we only need to solve for j=2 to nr-1
    % since the first and last points are fixed by the boundary condition
    sigma   =  D*dt/(dy(4)^2);
    for j   =   2:ny(4)-1
        Cnp14(j) =   sigma*Cn4(j+1)+(1-2*sigma)*Cn4(j)+sigma*Cn4(j-1);
    end
        Cnp14(ny(4)) =   Cn4(ny(4))+sigma*(Cn4(ny(4)-1)-Cn4(ny(4)));
    
        % 5 % Numerical solution - note that we only need to solve for j=2 to nr-1
    % since the first and last points are fixed by the boundary condition
    sigma   =  D*dt/(dy(5)^2);
    for j   =   2:ny(5)-1
        Cnp15(j) =   sigma*Cn5(j+1)+(1-2*sigma)*Cn5(j)+sigma*Cn5(j-1);
    end
        Cnp15(ny(5)) =   Cn5(ny(5))+sigma*(Cn5(ny(5)-1)-Cn5(ny(5)));
    
        % 6 % Numerical solution - note that we only need to solve for j=2 to nr-1
    % since the first and last points are fixed by the boundary condition
    sigma   =  D*dt/(dy(6)^2);
    for j   =   2:ny(6)-1
        Cnp16(j) =   sigma*Cn6(j+1)+(1-2*sigma)*Cn6(j)+sigma*Cn6(j-1);
    end
        Cnp16(ny(6)) =   Cn6(ny(6))+sigma*(Cn6(ny(6)-1)-Cn6(ny(6)));
    
    %Plot analytical vs numerical solution
    if plotFlag==1
        plot(y1,Cnp11, 'r+', y2,Cnp12, 'g+', y3,Cnp13, 'b+', y4,Cnp14, 'c+',  y5,Cnp15, 'm+', y6,Cnp16,'y+',  ya, CExact, 'k-')  
        xlabel('y (m)')
        ylabel('C (g/L)')
        ylim([0 50])
        title('Analytical vs Numerical Solution')
        h=legend('ny = 9','ny = 17','ny = 33','ny = 65','ny = 129','ny = 257','analytical');
        set(h,'location','northeast');
        pause(0.01);
    end
    
    % Copy solution to initial conditions for next iteration. Note only
    % copy the computed values otherwise you may overwrite the boundary
    % conditions!
    Cn1=Cnp11;  
    Cn2=Cnp12;
    Cn3=Cnp13;
    Cn4=Cnp14;
    Cn5=Cnp15;
    Cn6=Cnp16;
    
end


