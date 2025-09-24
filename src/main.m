% initialise model
init;

% main time stepping loop
while time-dt <= tend*yr
    
    % store previous position
    Xo = X;
    
    % calculate new orbital velocity
    for nj = 1:N
        D  = sum((X-X(nj,:)).^2,2).^0.5 + 1e-16;
        Fj = - (M.*M(nj))./D.^2 .* (X-X(nj,:))./D;
        V  = V + Fj./M .* dt;
    end
    V = V - V(1,:);
    
    X  = Xo + V.*dt;                         % update position of all bodies

    % detect collisions
    nj = 1;
    while nj < N
        Do      = sum((Xo-Xo(nj,:)).^2,2).^0.5 + 1e-16;
        D       = sum((X -X (nj,:)).^2,2).^0.5 + 1e-16;
        ind     = (D+Do)./2 < cls;
        ind(nj) = 0;
        ind     = find(ind>0);
        
        % merge collided bodies
        for i = 1:length(ind)
            X(nj,:)   = ((M(ind(i)).*X(ind(i),:)) + M(nj).*X(nj,:))./(M(ind(i))+M(nj));
            V(nj,:)   = ((M(ind(i)).*V(ind(i),:)) + M(nj).*V(nj,:))./(M(ind(i))+M(nj));
            C(nj,:)   = ((M(ind(i)).*C(ind(i),:)) + M(nj).*C(nj,:))./(M(ind(i))+M(nj));
            M(nj)     = M(ind(i))+ M(nj);
        end
        
        % remove consumed bodies
        X(ind,:)  = [];
        Xo(ind,:) = [];
        V(ind,:)  = [];
        C(ind,:)  = [];
        M(ind)    = [];
        N         = N-length(ind);

        CLS = CLS + length(ind); % update collision counter
        nj  = nj+1; % move to check collisions for next body
    end
    
    % update radial distance to sun
    r  = sum((X-X(1,:)).^2,2).^0.5 + 1e-32;  

    % update radii
    Rtot = sum(M.*C(:,1:3)./[2,1,0.5]+1e-16,2).^(1/3);
    Rrck = sum(M.*C(:,1:2)./[2,1    ]+1e-16,2).^(1/3);
    Rmtl = sum(M.*C(:,1:1)./[2      ]+1e-16,2).^(1/3);
    Rice = Rtot - Rrck - Rmtl;

    if ~mod(k,10); fprintf(1,'   -- %d;  time = %4.2f years;  bodies = %d;  collisions = %d\n',k,time/yr,N,CLS); end
    
    % plot and print model progress
    if ~mod(k,nop) || time-dt >= tend*yr
        output;
    end
    
    % update time and step count
    time = time + dt;
    k    = k + 1;
    
end