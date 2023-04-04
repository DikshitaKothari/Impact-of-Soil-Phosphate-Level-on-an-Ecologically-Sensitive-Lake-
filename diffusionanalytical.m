clc 
clear all 
close all 

timeend = 1000000;    % Ending time
nt = 101;             % Number of temoral points
plotFlag = 1;         % If want to plot

% Spatial (y) parameters
L      = 80000;                     % Size of the domain
ystart = 0;                         % Start of computational domain (m)
yend   = L;                         % End of computational domain (m)
ny     = 101;                       % Number of grid points along the  vector
dy     = (yend-ystart)/(ny-1);      % Spatial step size, delta y (m)
y      = ystart:dy:yend;            % Vector of grid points in y

% Temporal (t) parameters 
timestart = 0;                              % Starting time (0 seconds)
dt        = (timeend-timestart)/(nt-1);     % Temporal step size, delta t (s)
time      = timestart:dt:timeend;           % Vector of times

% Phyiscal parameters
phi0     = 50;                          % Temperature at y=0
phi_1    = 0;                           % Temperature at y=L
D        = 2.37*10^2                    % Diffusion equation constant                 
d        = sqrt(D)                      % Diffusion equation constant
toxicity = 0;                           % Toxic flag

% Calculate B coefficients using a for loop
nfs = 1000;               % Number of fourier terms               
B = zeros(1,nfs);         % Initialise B vector   
lambda=zeros(1,nfs);      % Initialise lambda vector
p = zeros(1,nfs);         % Initialise p vector
for n = 1:nfs;            % Compute values
    p(n) = (n-0.5)*pi/L;
    B(n)= -(2*phi0/(L*p(n)));                   
    lambda(n)=d*p(n);
end


%% Solve for C using three for loops in time, space, and Fourier series

% Loop through time
for i = 1:nt               % For each time
    
    t = time(i);           % time in seconds
    
    % Vector of zeros to initialise the Fourier series solution.
    % This should be re-initialised at each new time step.
    C = zeros(1,ny); 
    
    % Loop through space
    for j = 1:ny         % For each grid point        
        C(j) = phi0;         % Add the steady state solution
        
        % Loop through the Fourier series
        for n = 1:nfs
            C(j)  = C(j) + B(n)*sin(p(n)*y(j))*exp(-lambda(n)^2*t);  % Calculate series sum for T at y(i) and t
        end
    end
    
    % Plot the solution if requested
    if plotFlag==1
        plot(y,C); 
        ylim([0 51])
        grid on;                 
        xlabel('y (m)');            % Do not forget the unit from the plot
        ylabel('C (g/L)');
        title('Diffusion in the river')
        pause(0.001);               % Animate the plot (use instead of pause)'
    end

    if toxicity == 0            % Check toxicity flag 
        if C(ny) <= 50*10^-3    % Check against limit
            toxicity = 0
        else 
            toxicity = 1        % River is toxic
            T = t;              % Save time
        end
    end
end





