clear all;
close all;


%% parameters
% Coriolis/rotation parameters:
Omega = 7.3E-5; % angular rotation rate of Earth
latitude = 45; % degrees
f0 = 2*Omega*sind(latitude); % local Coriolis frequency

% geostrophic flow parameters:
dugdy = .75.*f0; % Ro_g = -dugdy/f0
F = sqrt(f0.*(f0-dugdy)); % effective Coriolis frequency

% wind stress parameters:
% [see eqns (15)-(16) and Fig. 1 of Whitt and Thomas 2015]
taumaga = .06; % N/m^2
taumagb = 0;   % N/m^2
thetaw = 0;
rho = 1.025e3; % reference density of sea water kg/m^3
omegaw = .5.*f0; % wind forcing angular frequency, s^{-1} 
% defined as in Whitt and Thomas 2015

% MLD is constant here:
H = 25; % meters

% Slab mixed layer damping parameter:
r = 1./(2.*24.*3600); % s^{-1}

% initial conditions for analytic solution (see text following
% (B10) and (B19) in WT15:
% change these two:
u0 = 0;
v0 = 0;
% but, do not change these two:
dvdt0 = -r.*v0 - f0.*u0 -taumaga.*sin(thetaw)./(rho.*H);
dudt0 = -r.*u0 + (F.^2)./f0.*v0 + taumaga.*cos(thetaw)./(rho.*H);

% find numerical solution to (1)-(2) in WT15 for comparison with the
% analytic solution: 
timepts = linspace(0,15.*2.*pi./F,10000);
uvec = [u0 v0]';
options = odeset('RelTol',1e-4,'AbsTol',[1e-5 1e-5]);
[tout,uvout] = ode45(@(tout,uvout) odefunhelperTEab(tout,uvout,omegaw,taumaga,taumagb,thetaw,H,r,rho,F,f0),timepts,uvec,options);
t2 = tout(:,1)./(2.*pi./F); 
%% find the analytic solution under sustained forcing (following
% appendix B of WT15) :

% first find the particular/forced solution:

X = F.^2 + r.^2 - omegaw.^2; % (B6) 
Y = 2.*r.*omegaw;  % Y in WT15; (B7) 
SV = taumaga.*omegaw.*sin(thetaw) - taumagb.*r.*cos(thetaw) + taumagb.*f0.*sin(thetaw); % (B5)
CV = -(taumagb.*omegaw.*cos(thetaw) + r.*taumaga.*sin(thetaw) + f0.*taumaga.*cos(thetaw)); %  (B4) 

SU = -taumaga.*omegaw.*cos(thetaw) - r.*taumagb.*sin(thetaw) - ((F.^2)./f0).*taumagb.*cos(thetaw);  % (B16) 
CU = -taumagb.*omegaw.*sin(thetaw) + r.*taumaga.*cos(thetaw) - ((F.^2)./f0).*taumaga.*sin(thetaw); % (B15) 

phinow = angle((CV.*X - Y.*SV).*1i +(SV.*X + CV.*Y)); % (B3)

phinow2 = angle((CU.*X - Y.*SU).*1i +(SU.*X + CU.*Y)); % (B14)

AV = (1./(rho.*H)).*sqrt(SV.^2 + CV.^2)./(sqrt((X).^2 + (Y).^2)); % (B2)
AU =  (1./(rho.*H)).*sqrt(SU.^2 + CU.^2)./(sqrt((X).^2 + (Y).^2)); % (B13)

% then find the transient/unforced solution:
C2 = -1./(F).*(f0.*u0 + AV.*(omegaw.*cos(phinow) + r.*sin(phinow)) + taumaga.*sin(thetaw)./(rho.*H)); % (B10) 
C1 = v0-AV.*sin(phinow); % (B9) 

D2 = -1./(F).*(-F.^2./f0.*v0 + AU.*(omegaw.*cos(phinow2) + r.*sin(phinow2)) - taumaga.*cos(thetaw)./(rho.*H));  % (B19) 
D1 = u0-AU.*sin(phinow2); % (B18) 

% then take the sum to obtain the full solution = transient/homogeneous + particular/forced 
Utota = AU.*sin(omegaw.*tout(:,1) + phinow2) + D2.*exp(-r.*tout(:,1)).*sin(F.*tout(:,1)) + D1.*exp(-r.*tout(:,1)).*cos(F.*tout(:,1));
Vtota = AV.*sin(omegaw.*tout(:,1) + phinow) + C2.*exp(-r.*tout(:,1)).*sin(F.*tout(:,1)) + C1.*exp(-r.*tout(:,1)).*cos(F.*tout(:,1));


%% calculating the response function (18) [see (B20)-(B23)]:
deltatheta = phinow - phinow2;
% (B21):
lengthone = sqrt(1./((AU.^2 - (AU.^4 + AV.^4 - 2.*AU.^2.*AV.^2 + 4.*AU.^2.*AV.^2.*cos(deltatheta).^2).^(1./2) + AV.^2)./(2.*AU.^2.*AV.^2.*sin(deltatheta).^2)));
lengthtwo = sqrt(1./((AU.^2 + (AU.^4 + AV.^4 - 2.*AU.^2.*AV.^2 + 4.*AU.^2.*AV.^2.*cos(deltatheta).^2).^(1./2) + AV.^2)./(2.*AU.^2.*AV.^2.*sin(deltatheta).^2)));
La = max(lengthone,lengthtwo);
Lb = min(lengthone,lengthtwo);
Ta = max(taumaga,taumagb);
Tb = min(taumaga,taumagb);
display('Testing response fn')
% compute the averaged velocity response function (B23), which
% appears in the numerator of (18), numerically and then
% analytically (the two should be nearly identical):
testingresponse1 = nanmean(abs(AU.*sin(linspace(0,4.*pi,1e6) + phinow2) + 1i.*AV.*sin(linspace(0,4.*pi,1e6) + phinow)))
Ubra = (La./(2.*pi)).*(ellipticE(2.*pi,1-(Lb./La).^2))
% compute the averaged forcing function (17), 
% which appears in the denominator of (18) 
% numerically and analytically (the two should be nearly identical):
testingresponse2 = nanmean(abs(Ta.*sin(linspace(0,4.*pi,1e6)) + 1i.*Tb.*cos(linspace(0,4.*pi,1e6))))
Taubra = (Ta./(2.*pi)).*(ellipticE(2.*pi,1-(Tb./Ta).^2))

% Eqn (18):
responsefn = rho.*H.*f0.*Ubra./Taubra
%% plot numerical and analytical solutions for the
% sustained forcing case (they should be the same)
figure
subplot(1,2,1),...
plot(t2,uvout(:,1),t2,uvout(:,2),t2,Utota,'--',t2,Vtota,'--','linewidth',2);
xlim([0 15])
ylim([-1 1])
xlabel('Effective Inertial Periods (Ft/2\pi)')
ylabel('[m/s]')
grid on
legend('u_n','v_n','u_a','v_a')
title('Sustained forcing: U_{ML},V_{ML}')
% calculate energy budget (19) for the sustained forcing case:

% wind stress (15)-(16):
tauxnow = taumaga.*cos(omegaw.*tout(:,1)).*cos(thetaw) - taumagb.*sin(omegaw.*tout(:,1)).*sin(thetaw);
tauynow = -taumaga.*cos(omegaw.*tout(:,1)).*sin(thetaw) - taumagb.*sin(omegaw.*tout(:,1)).*cos(thetaw);

% left hand side (LHS) of (19)
E = .5.*(Utota.^2 + Vtota.^2);
dt = tout(2,1)-tout(1,1);
% LSP on RHS of (19):
cross = cumsum(dt.*(-dugdy.*Utota.*Vtota),1);
% DAMP on RHS of (19):
damp = cumsum(dt.*(-r.*Utota.^2 - r.*Vtota.^2),1);
% WORK on RHS of (19):
force = cumsum(dt.*(tauxnow.*Utota./(rho.*H) + tauynow.*Vtota./(rho.*H)),1);
% sum of all three terms on the right hand side (RHS) of (19):
rhs = cumsum(dt.*(-dugdy.*Utota.*Vtota - r.*Utota.^2 - r.*Vtota.^2 + tauxnow.*Utota./(rho.*H) + tauynow.*Vtota./(rho.*H)),1);

subplot(1,2,2),...
plot(t2,E,t2,rhs,'--',t2,damp,t2,force,t2,cross,'--','linewidth',2);
legend('E','rhs','damp','force','cross')
xlabel('Effective Inertial Periods (Ft/2\pi)')
ylabel('[m^2/s^2]')
grid on
xlim([0 15])
ylim([-1 1])
title('Sustained forcing: E_{ML} budget')


%% Calculate the viscous decay scenario as in Fig. 7
% after ndays of forcing
ndays = 1; % 1 day of forcing as in Fig. 7 
[C,stidx] = min(abs(tout - ndays.*3600.*24));

U02 = Utota(stidx);
V02 = Vtota(stidx);
DUDT02 = ((F.^2)./f0).*V02 - r.*U02;
DVDT02 = -f0.*U02 - r.*V02;
V3 = zeros(size(Vtota));
U3 = zeros(size(Vtota));
taux2 = tauxnow;
taux2(stidx:end) = 0;
tauy2 = tauynow;
tauy2(stidx:end) = 0;

V3(1:stidx-1) = Vtota(1:stidx-1);
V3(stidx:end) = V02.*exp(-r.*(tout(stidx:end,1) - tout(stidx,1))).*cos(F.*(tout(stidx:end,1) - tout(stidx,1))) + ...
    ((DVDT02 + r.*V02)./F).*exp(-r.*(tout(stidx:end,1) - tout(stidx,1))).*sin(F.*(tout(stidx:end,1) - tout(stidx,1)));
U3(1:stidx-1) = Utota(1:stidx-1);
U3(stidx:end) = U02.*exp(-r.*(tout(stidx:end,1) - tout(stidx,1))).*cos(F.*(tout(stidx:end,1) - tout(stidx,1))) + ...
    ((DUDT02 + r.*U02)./F).*exp(-r.*(tout(stidx:end,1) - tout(stidx,1))).*sin(F.*(tout(stidx:end,1) - tout(stidx,1)));

figure;
subplot(1,2,1),...
plot(tout./86400,U3,tout./86400,V3);
xlabel('[Days]')
ylabel('[m/s]')
legend('U_{ML}','V_{ML}');
xlim([0 5])
grid on
title('Momentum in the viscous decay')

% time-integrated energy budget:
E2 = .5.*(U3.^2 + V3.^2);
dt = tout(2,1)-tout(1,1);
cross2 = cumsum(dt.*(-dugdy.*U3.*V3),1);
damp2 = cumsum(dt.*(-r.*U3.^2 - r.*V3.^2),1);
force2 = cumsum(dt.*(taux2.*U3./(rho.*H) + tauy2.*V3./(rho.*H)),1);
rhs2 = cumsum(dt.*(-dugdy.*U3.*V3 - r.*U3.^2 - r.*V3.^2 + taux2.*U3./(rho.*H) + tauy2.*V3./(rho.*H)),1);

subplot(1,2,2),...
plot(tout./86400,E2,tout./86400,rhs2,'--',t2,damp2,t2,force2,t2,cross2,'--','linewidth',2);
legend('E','rhs','damp','force','cross')
xlabel('[Days]')
xlim([0 5])
ylabel('[m^2/s^2]')
title('Time-integrated energy budget for viscous decay (Fig 7a)')
grid on


