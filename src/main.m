% initialise model
init;


% main time stepping loop
while time <= tend*yr

    %--------------------------------------------------------------
    % Store previous state
    %--------------------------------------------------------------
    Xo = X;
    Vo = V;
    Fjo = Fj;

    %--------------------------------------------------------------
    % Velocity Verlet
    %--------------------------------------------------------------
    
    % update positions
    Ao = Fjo ./ M;

    X = Xo + Vo*dt + 0.5*Ao*dt^2;

    % Gravitational force at new positions
    [Fj] = get_forces(M,X,N);

    % update velocities 
    V = Vo + 0.5*(Fjo + Fj)./M * dt;

    %--------------------------------------------------------------
    % Collision detection and merging
    %--------------------------------------------------------------
    nj = 1;

    while nj < N  % loop through all bodies

        % Relative position at beginning of timestep
        r0 = Xo - Xo(nj,:);

        % Change in relative position over timestep
        dr = (X - X(nj,:)) - r0;

        % Closest approach assuming linear trajectories
        tau = - sum(r0.*dr,2) ./ (sum(dr.^2,2)+eps);
        tau = max(0,min(1,tau));

        % Minimum separation
        rmin = r0 + tau.*dr;
        dmin = sqrt(sum(rmin.^2,2));

        % Detect collisions
        ind = dmin < cls;
        ind(nj) = false;

        if any(ind)

            % Conservative merger
            Mnew    = M(nj) + sum(M(ind));
            X(nj,:) = (M(nj)*X(nj,:) + sum(M(ind).*X(ind,:),1)) / Mnew;
            V(nj,:) = (M(nj)*V(nj,:) + sum(M(ind).*V(ind,:),1)) / Mnew;
            C(nj,:) = (M(nj)*C(nj,:) + sum(M(ind).*C(ind,:),1)) / Mnew;
            M(nj)   = Mnew;

            % Remove consumed bodies
            X(ind,:)   = [];
            V(ind,:)   = [];
            Fj(ind,:)  = [];
            Xo(ind,:)  = [];
            Vo(ind,:)  = [];
            Fjo(ind,:) = [];
            C(ind,:)   = [];
            M(ind)     = [];

            % update number of bodies post collision
            N = N - sum(ind);

            % update number of collisions
            CLS = CLS + sum(ind);
        end

        nj = nj + 1;
    end

    %--------------------------------------------------------------
    % Update model state
    %--------------------------------------------------------------

    % Keep Sun stationary
    X = X - X(1,:);

    % Update radial distance to Sun
    r = sqrt(sum((X-X(1,:)).^2,2)) + eps^2;

    % Update radii
    [Rtot,Rmtl,Rrck,Rsun,Rggt] = get_radii(M,C);

    %--------------------------------------------------------------
    % Report Model progress
    %--------------------------------------------------------------
    if ~mod(step,10)
        fprintf(1,...
            '   -- %d;  time = %4.2f years;  bodies = %d;  collisions = %d\n',...
            step,time/yr,N,CLS);
    end

    % produce output figures, save model state
    if ~mod(step,nop) || time >= tend*yr
        output;
    end

    % increment time and step count
    time = time + dt;
    step = step + 1;

end