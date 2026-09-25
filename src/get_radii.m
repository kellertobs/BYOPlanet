function [Rtot,Rmtl,Rrck,Rsun,Rggt] = get_radii(M,C)

Rtot = sum(M.*C(:,1:3)./[1.6,0.8,0.4]+eps,2).^(1/3);
Rrck = sum(M.*C(:,1:2)./[1.6,0.8    ]+eps,2).^(1/3);
Rmtl = sum(M.*C(:,1:1)./[1.6        ]+eps,2).^(1/3);
Rsun = Rtot(1);
Rggt = Rtot(2);

end