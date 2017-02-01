function [dudt] = odefunhelperTEab(t,u,sigma,taumaga,taumagb,thetaw,Z,c,rho,F,f0)
dudt = zeros(2,1);
tauxnow = taumaga.*cos(sigma.*t).*cos(thetaw) - taumagb.*sin(sigma.*t).*sin(thetaw);
tauynow = -taumaga.*cos(sigma.*t).*sin(thetaw) - taumagb.*sin(sigma.*t).*cos(thetaw);
dudt(1) =  ((F.^2)./f0).*u(2) - c.*u(1) + (tauxnow./(rho.*Z));
dudt(2) = -f0.*u(1) - c.*u(2) + (tauynow./(rho.*Z));
end