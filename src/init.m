% create output directory
if ~exist(['../out/',runID],'dir'); mkdir(['../out/',runID]); end

if restart

if     restart < 0  % restart from last continuation frame
    name = ['../out/',runID,'/',runID,'_cont.mat'];
elseif restart > 0  % restart from specified continuation frame
    name = ['../out/',runID,'/',runID,'_',num2str(restart),'.mat'];
end
if exist(name,'file')
    fprintf('\n  restart from %s \n\n',name);
    load(name,'N','M','C','X','V','Fj','CLS','time','step');

    % update radii
    Rear = 1;
    [Rtot,Rmtl,Rrck,Rsun,Rggt] = get_radii(M,C);


    % calculate duration of 1 year (orbital period for body of M = 1, R = 1)
    yr = 2*pi/sqrt(M(1));

    % calculate time step size
    dt = yr/100;

    % update time and step count
    time = time + dt;
    step = step + 1;

else
    restart = 0; 
end

end

if ~restart

% initialise random number generator
rng(seed);

% initialise body mass and position
M  = [MStr; MGgt; min(MPls*100,max(MPls/100,normrnd(MPls,MPls/2,N-2,1)))];  % body mass
X  = [0,0,0; 5,0,0; randn(N-2,3).*[7,7,0]];  % body position
r  = sum((X-X(1,:)).^2,2).^0.5 + eps^2;  % radial distance to sun
X(3:end,3) = sind(Incl).*r(3:end).*randn(N-2,1);

[Fj] = get_forces(M,X,N);

% initialise body composition (metal = 1, rock = 2; ice = 3)
Cmtl = min(1,max(0,0.65-r.^0.25/3.5) .* (1 + randn(N,1)./100));
Crck = 1-Cmtl;
Cice = min(10,max(0,max(0,r-5).^0.75) .* (1 + randn(N,1)./100));
C    = [Cmtl,Crck,Cice]./(Cmtl+Crck+Cice);
C(1,:) = [0 0 1];
C(2,:) = [0.1 0.2 0.7];

% calculate radii
Rear = 1;  % Earth radius
[Rtot,Rmtl,Rrck,Rsun,Rggt] = get_radii(M,C);

% calculate initial orbital velocity
R  = (X-X(1,:))./r;
V  = sqrt(M(1)./r) .* [R(:,2),-R(:,1),R(:,3)] .* (1+[0 0 0;0 0 0;randn(N-2,3)./10]);

% calculate duration of 1 year (orbital period for body of M = 1, R = 1)
yr = 2*pi/sqrt(M(1));

% calculate time step size
dt = yr/100;

CLS  = 0;  % initialise collision count
time = 0;  % initialise time
step = 0;  % initialise time step counter

end

Vo   = V;
Xo   = X;

% print start of simulation
fprintf(1,'\n\n*****  Start B-Y-O Planet Simulation  *****\n\n')