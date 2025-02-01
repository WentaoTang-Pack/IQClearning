function v = twophase_model(States, Inputs, Disturbances)
% Parameters and nominal values
cp = 0.08; rho = 15000;
EaR = 13230.695; k0 = 1.0e11;
TA0 = 310; TB0_nom = 298; FA0 = 171.25;
FB0_nom = 300; FL_nom = 375;
Vtot = 3; DHR = -50; DHV = 20;
FV_nom = 50; Q_nom = -3422.6;
% Saturated pressures
% PAsat = @(T) exp(30.5-3919.7/(T-34.1));
% PCsat = @(T) exp(30.0-5000.0/(T+70.0)); 
dlnPAsat = @(T) 3919.7/(T-34.1)^2; 
dlnPCsat = @(T) 5000.0/(T+70.0)^2; 
% Steady state (u=0)
xss = [12834.252422604806; 0.7157887166210231; 12812.88206760394; ...
    0.2378948377838636; 0.6766666666666666; 341.50882352941176];
% States
MV = xss(1) + States(1); yA = xss(2) + States(2); ML = xss(3) + States(3); 
xA = xss(4) + States(4); xB = xss(5) + States(5); T = xss(6) + States(6); 
FV = FV_nom + Inputs(1)*1; FL = FL_nom - Inputs(2)*1; Q = Q_nom + Inputs(3)*25; 
FB0 = FB0_nom + Disturbances(1)*1; TB0 = TB0_nom + Disturbances(2)*2.5; 
rtot = k0*exp(-EaR/T)*ML*rho*xA*xB; 
dMV_0 = FA0 - FV; dMV_1 = -1; dMV_2 = 1; dyA_0 = FA0*(1-yA)/MV; 
dyA_1 = -(1-yA)/MV; dyA_2 = -yA/MV; dML_0 = FB0 - rtot - FL; dML_1 = 1; dML_2 = -1;
dxA_0 = -(rtot*(1-xA) + FB0*xA)/ML; dxA_1 = (1-xA)/ML; dxA_2 = xA/ML;
dxB_0 = (FB0 - rtot)*(1-xB)/ML; dxB_1 = -xB/ML; dxB_2 = xB/ML; 
dT_0 = 1/(ML+MV)/cp * (FA0*cp*(TA0-T)+FB0*cp*(TB0-T) - rtot*(DHR-cp*(T-298.0)) + Q);
dT_1 = 1/(ML+MV)/cp * DHV;
dT_2 = -1/(ML+MV)/cp * DHV; 
eta00 = dlnPAsat(T)*dT_0 + 1/xA*dxA_0 - 1/MV*dMV_0-1/T*dT_0-1/(rho*Vtot-ML)*dML_0 - 1/yA*dyA_0;
eta01 = dlnPAsat(T)*dT_1 + 1/xA*dxA_1 - 1/MV*dMV_1-1/T*dT_1-1/(rho*Vtot-ML)*dML_1 - 1/yA*dyA_1;
eta02 = dlnPAsat(T)*dT_2 + 1/xA*dxA_2 - 1/MV*dMV_2-1/T*dT_2-1/(rho*Vtot-ML)*dML_2 - 1/yA*dyA_2; 
eta10 = dlnPCsat(T)*dT_0 - 1/(1-xA-xB)*(dxA_0+dxB_0) - 1/MV*dMV_0-1/T*dT_0-1/(rho*Vtot-ML)*dML_0 + 1/(1-yA)*dyA_0; 
eta11 = dlnPCsat(T)*dT_1 - 1/(1-xA-xB)*(dxA_1+dxB_1) - 1/MV*dMV_1-1/T*dT_1-1/(rho*Vtot-ML)*dML_1 + 1/(1-yA)*dyA_1; 
eta12 = dlnPCsat(T)*dT_2 - 1/(1-xA-xB)*(dxA_2+dxB_2) - 1/MV*dMV_2-1/T*dT_2-1/(rho*Vtot-ML)*dML_2 + 1/(1-yA)*dyA_2; 
NA = -1/(eta01*eta12-eta02*eta11) * (eta12*eta00 - eta02*eta10); 
NC = -1/(eta01*eta12-eta02*eta11) * (eta01*eta10 - eta11*eta00); 
dMV = dMV_0 + dMV_1*NA + dMV_2*NC; 
dyA = dyA_0 + dyA_1*NA + dyA_2*NC; 
dML = dML_0 + dML_1*NA + dML_2*NC; 
dxA = dxA_0 + dxA_1*NA + dxA_2*NC; 
dxB = dxB_0 + dxB_1*NA + dxB_2*NC;
dT = dT_0 + dT_1*NA + dT_2*NC; 
v = [dMV; dyA; dML; dxA; dxB; dT]*60;
end

