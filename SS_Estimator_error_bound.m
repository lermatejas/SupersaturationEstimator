%% -------------INAUT-Naciional University of San Juan ------------------
% This algorithm simulates the supersaturation estimator with error bounds. 
% The authors of this study are:
% 1-Phd(c) Ing. Humberto Morales, 
% 2-Dr. Ing. Fernando Di Sciascio, 
% 3-Dra. Ing. Adriana Amicarrelli, 
% 4-Msc. Ing. Estefania Aguirre.
% The estimator is based on the Crystalization C model proposed by Damour 
% and the Crystalization model proposed by Lajos. 
% The authors have implemented an own methodology to obtain and model based estimator
% and it error bound of  3 standard deviations
% Version 1.0
%%------------------------------------------------------------------------
%% 
clear all; clc ;  warning off,
%% ---Times for simulations--------------------------------------
Ts=30;                 % Sample Time
TimeS=4800;            % Simulation Time 
global  ts2            % Time for Rate Transition blokof simulink
ts2=30;
ts=30;
%% ------ Load Model Prameters----------------------------------
p=param();             
%% -------Calculate the inicial conditios of state variables:-----
% Mw: mass water;Mi: mass impurities;Mc: mass of crystals; Ms: mass sucrose; (kg)
% T: temperature (C); 
% MTo:  Initial Total Mass (kg); dens: density (kg/m3); Vpie: initial Volume (m3)
% cc_o: Initial Mass of crystals (%);  T0: Initial Temperature (C)
%-----------------------------------------------------------------------------
MTo=p.dens*p.Vpie;   
syms Ms0 Mi0 Mw0
Mc0=p.cc_o*MTo;                     
Sol = solve([100*(Ms0+Mi0)/(Ms0+Mi0+Mw0)==p.Bx_ml, 100*(Ms0+Mi0+Mc0)/(Ms0+Mi0+...
    Mc0+Mw0)==p.Bx_mg,100*(Ms0+Mc0)/(Ms0+Mi0+Mc0)==p.P_mg], [Ms0,Mi0,Mw0]);
Ms0=double(Sol.Ms0);
Mi0=double(Sol.Mi0);
Mw0=double(Sol.Mw0);
T0=p.T0;                               
%-----------------------CI Filter--------------------------------------------
 Mw0f=Mw0-250;Mi0f=Mi0-250;Mc0f=Mc0-250;Ms0f=Ms0-250;T0f=T0;   
  CIf=[Mw0f;Mi0f;Mc0f;Ms0f;T0f];
  
%% --Block that generates a reference for simulation with the PID(2DOF) implemented in Simulink; 
%    to use set the switch  that selects between PID (2DOF) and Acc_ctrl in the simulink file --

%---Variable Reference-------------------------
% SS_ini(1,1)=1.26;
% for i=1:1:TimeS/Ts
%     SS_ini(1,i)=SS_ini(1,i)-0.015+0.03*rand;
%       if SS_ini(1,i)<1.258
%         SS_ini(1,i)=1.258;
%       end
%       if SS_ini(1,i)>1.262
%          SS_ini(1,i)=1.262;
%       end
%     Ref_ss(i)=SS_ini(1,i);
%     SS_ini(1,i+1)=Ref_ss(i);
% end

%--Fixed reference-----
Ref_ss=ones(161,1)*1.26;
 Ref_ss=Ref_ss';

%---Generate your own control actions (without using PID) to simulate; 
% to use set the switch  that selects between PID (2DOF) and Acc_ctrlin the simulink file --
for i=1:1:TimeS/Ts
Acc_ctrl(i)=0.0056-0.00001*2*i;
end
Acc_ctrl=Acc_ctrl'; %(m3/s)

%% ------Mass of evaporate vapor------------
% Use this typical value or provide your own data vector at each time step; 
% this value is associated with the Step block's Evaporated Vapor Mass in Simulink
M_vapor=1.34;  % masa vapor (kg/seg)

%% Simulation; load Simulink File 
load_system('Evaporator') %Invisibly load Simulink model
sim('Evaporator')
%-----Simulink returns the model states -and the supersaturation 
C=Pxe.Data*100;                                                     % Covariance Matrix
Mw = Mw.Data; Mi = Mi.Data; Mc = Mc.Data; Ms = Ms.Data; T = T.Data; % States 
x=[Mw Mi Mc Ms T];                                                  % State vectors
time=SS.Time/3600;                                                  % Time Vector
%% This block calculates the sensitivity functions (partial derivatives of SS with respect to the states
syms x1 x2 x3 x4 x5
SS_sym=(x4/x1)/(((64.447+0.08222*x5+1.66169e-3*x5^2-1.558e-6*x5^3-...
    4.63e-8*x5^4)*(-0.067*x2/x1+0.96+0.04*exp(-2.8*x2/x1))/...
    (100-(64.447+0.08222*x5+1.66169e-3*x5^2-1.558e-6*x5^3-4.63e-8*x5^4))));
DSSDx1 = diff(SS_sym,x1); DSSDx2 = diff(SS_sym,x2); DSSDx3 = diff(SS_sym,x3);
DSSDx4 = diff(SS_sym,x4); DSSDx5 = diff(SS_sym,x5);

dSSdx1_sym=[];dSSdx2_sym=[];dSSdx3_sym=[];dSSdx4_sym=[];dSSdx5_sym=[];
for k=1:161
% Symbolic partial derivatives.
dSSdx1_sym(k)=subs(subs(subs(subs(subs(DSSDx1,x1,Mw(k)),x2,Mi(k)),x3,Mc(k)),x4,Ms(k)),x5,T(k));
dSSdx2_sym(k)=subs(subs(subs(subs(subs(DSSDx2,x1,Mw(k)),x2,Mi(k)),x3,Mc(k)),x4,Ms(k)),x5,T(k));
dSSdx3_sym(k)=subs(subs(subs(subs(subs(DSSDx3,x1,Mw(k)),x2,Mi(k)),x3,Mc(k)),x4,Ms(k)),x5,T(k)); % 0 ya que SS no depende de x3 (Mc)
dSSdx4_sym(k)=subs(subs(subs(subs(subs(DSSDx4,x1,Mw(k)),x2,Mi(k)),x3,Mc(k)),x4,Ms(k)),x5,T(k));
dSSdx5_sym(k)=subs(subs(subs(subs(subs(DSSDx5,x1,Mw(k)),x2,Mi(k)),x3,Mc(k)),x4,Ms(k)),x5,T(k));
% Numerical partial derivatives.
dSSdMw(k)=double(dSSdx1_sym(k));dSSdMi(k)=double(dSSdx2_sym(k));dSSdMc(k)=double(dSSdx3_sym(k)); 
dSSdMs(k)=double(dSSdx4_sym(k)); dSSdT(k)=double(dSSdx5_sym(k));
end

% Theoretically, the partial derivatives should be constant because SS is not an explicit function of time. 
% However, in simulations, they may vary slightly over time. We can consider them as constant
dSSdMw=mean(dSSdMw); dSSdMi=mean(dSSdMi); dSSdMc=mean(dSSdMc);
dSSdMs=mean(dSSdMs);dSSdT=mean(dSSdT);
% We used this values
dSSdMw=-2.610e-04; dSSdMi=1.83e-05; dSSdMc=0; dSSdMs=7.04e-05; dSSdT=-0.0152;

%% -----First order Taylor approximation.------------------------------------
xe=xe.Data; % Vector of estimated states from Simulink.
for k=1:161
SSh(k)=(Ms(k)/Mw(k))/(((64.447+0.08222*T(k)+1.66169e-3*T(k)^2-1.558e-6*T(k)^3-...
    4.63e-8*T(k)^4)*(-0.067*Mi(k)/Mw(k)+0.96+0.04*exp(-2.8*Mi(k)/Mw(k)))/...
    (100-(64.447+0.08222*T(k)+1.66169e-3*T(k)^2-1.558e-6*T(k)^3-4.63e-8*T(k)^4))));
w1(k)=xe(k,1)-x(k,1); w2(k)=xe(k,2)-x(k,2); w4(k)=xe(k,4)-x(k,4); w5(k)=xe(k,5)-x(k,5);
 % SSe_Taylor(k) = SSh(k) + dSSdMw*w1(k) + dSSdMi*w2(k) + dSSdMs*w4(k) + dSSdT*w5(k);
   SSe_Taylor(k) = SSh(k)' - 2.608e-04*w1(k) + 1.8284e-05*w2(k) + 7.0386e-05*w4(k) - 0.0152*w5(k);
end
%% Variance and standar desviation calculation
varMw=[]; varMi=[]; varMc=[]; varMs=[]; VARt=[];
covarMwMi=[];covarMwMc=[];covarMwMs=[];covarMwT=[];covarMiMc=[];covarMiMs=[];covarMiT=[];covarMcMs=[];covarMcT=[];covarMsT=[];

for k=1:161
varMw(k)=C(1,1,k); varMi(k)=C(2,2,k); varMc(k)=C(3,3,k); varMs(k)=C(4,4,k);varT(k)=C(5,5,k);
covarMwMi(k)=C(1,2,k); covarMwMc(k)=C(1,3,k); covarMwMs(k)=C(1,4,k);
covarMwT(k)=C(1,5,k); covarMiMc(k)=C(2,3,k); covarMiMs(k)=C(2,4,k);
covarMiT(k)=C(2,5,k);covarMcMs(k)=C(3,4,k); covarMcT(k)=C(3,5,k);
covarMsT(k)=C(4,5,k);

varSS(k) = dSSdMw^2*varMw(k)+dSSdMi^2*varMi(k)+dSSdMs^2*varMs(k)+dSSdT^2*varT(k)+...
    2*(dSSdMw*dSSdMi*covarMwMi(k)+dSSdMw*dSSdMs*covarMwMs(k)+dSSdMw*dSSdT*...
    covarMwT(k)+dSSdMi*dSSdMs*covarMiMs(k)+dSSdMi*dSSdT*covarMiT(k)+dSSdMs*dSSdT*covarMsT(k));
varSS_solo_de_T(k) = dSSdT^2*varT(k);
end

stdSS=sqrt(varSS); 
% Three standard deviations include all the numbers for 99.7% of the sample population being studied
SS_up  = SSh + 3*stdSS;
SS_low = SSh - 3*stdSS;
%% Estimated Supersaturation from simulink
 SSe = SSe_con_T.Data;
 MSE= immse(SSh,SSe') % Mean Squared Error
%% Figure setting
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultlineMarkerSize',12)
set(0,'DefaultlineMarkerFace','b')
set(0,'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight','demi')
set(0,'DefaultTextFontSize', 20, 'DefaultTextFontWeight','demi')

%------ Figures------------------
figure(1)
X=[time' time(end:-1:1)']; Y=[SS_up SS_low(end:-1:1)];fill(X,Y,[0.85 0.85 0.85]);hold on; 
plot(time,SSh','-.k',time,SSe,'r',time,SSe_Taylor,'--b');hold on;plot(time,Ref_ss','-.c')
xlim([0 1.3]);grid on; legend('3\sigma_{SS} error band','SS model','Full SS estimator model',...
    'Taylor first order estimator','Reference','Location','best');
xlabel('Time(h)');ylabel('Supersaturation');title('Supersaturation estimator')

