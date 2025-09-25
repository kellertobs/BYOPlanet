% create output directory
if ~exist(['../out/',runID],'dir'); mkdir(['../out/',runID]); end

% initialise random number generator
rng(seed);

% initialise body mass and position
M  = [MStr; MGgt; min(MPls*100,max(MPls/100,normrnd(MPls,MPls/5,N-2,1)))];  % body mass
X  = [0,0,0; 5,0,0; randn(N-2,3).*[6,6,1/3]];  % body position
D  = squareform(pdist(X,'euclidean')) + eps;  % mutual distance matrix
r  = sum((X-X(1,:)).^2,2).^0.5 + eps^2;  % radial distance to sun

% initialise body composition (metal = 1, rock = 2; ice = 3)
Cmtl = min(1,max(0,0.5-r.^0.5./6 + randn(N,1)./100));
Crck = 1-Cmtl;
Cice = min(30,max(0,3.*max(0,r-5).^0.5) .* (1 + randn(N,1)./20));
C    = [Cmtl,Crck,Cice]./(Cmtl+Crck+Cice);

% calculate radii
Rtot = sum(M.*C(:,1:3)./[2,1,0.5]+eps,2).^(1/3);
Rrck = sum(M.*C(:,1:2)./[2,1    ]+eps,2).^(1/3);
Rmtl = sum(M.*C(:,1:1)./[2      ]+eps,2).^(1/3);
Rsun = (M(1)/0.5).^(1/3);
Rggt = (M(2)/0.5).^(1/3);
Rear =    (1/1.5).^(1/3);

% calculate initial orbital velocity
R  = (X-X(1,:))./r;
V  = sqrt(M(1)./r) .* [R(:,2),-R(:,1),R(:,3)] .* (1+randn(N,3).*[0 0 0;0 0 0;ones(N-2,1).*[0.1,0.1,0.5]]);

% calculate duration of 1 year (orbital period for body of M = 1, R = 1)
yr = 2*pi/sqrt(M(1));

% calculate time step size
dt = yr/100;

CLS  = 0;  % initialise collision count
time = 0;  % initialise time
k    = 0;  % initialise time step counter

% print start of simulation
fprintf(1,'\n\n*****  Start B-Y-O Planet Simulation  *****\n\n')