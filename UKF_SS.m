% Unscented Kalman Filter from Simon(2006)

function ze = UKF_SS(u)
Bx_ml_medida= u(1);
V_medida=u(2);
T_medida=u(3);

y=[Bx_ml_medida;V_medida;T_medida];  % Measurement

Flujo=u(4);         % Control Action
Mvapor_f=u(5);      % Control Action
xee = u(6:10);      %  UFK States
%%
%Simbolicas
persistent Qe Re PP n
global aux2 ts2
% syms Mw_f Mi_f Mc_f Ms_f T_f
 ts2=30;
%---------------UKF Initialization-------------------------------
if isempty(aux2)
    v1 = [4 2 4.5 2.5 1];           % Weights of the main diagonal in PP
    v2 = [1 1 1 1 1];               % Weights of the main diagonal in Qe
    v3 = [1 1 1 ];

    PP = diag(v1);
    Qe = 0.01*diag(v2);
    Re = 0.0001*diag(v3);

    aux2 = 1;
    n = 5; % System Order
end
%% Model parameters
p=param(); %load 
%% *************************************************************************
% UKF ALGORITHM Equations466 4.56 to 14.67 of Simon [54].
     %Generation of sigma points: 1st set
%**************************************************************************
% Sigma Points as in Wan, E. A. and Van der Merwe, R. The Unscented Kalman
% Filter for Nonlinear Estimation
%--------------------------(3a)-------------------------------------
%Sigma Point Selections 
% sqPP = square root of (n*PP). The Cholesky factorization is used.
sqPP =  (chol(n*PP))';

for i=1:2*n
    if i<=n
        xhat1(:,i) = xee + sqPP(:,i);
    else
        xhat1(:,i) = xee - sqPP(:,i-n);
    end
end

% Sample time (discretization)
ts = ts2; % [min]
%---------------------(3b)-------------------------------------------------
% Sigma points transformation through the nonlinear map
% Se usa la ecuaci?n no lineal del  sistema para  transformar los 
% sigma-puntos en los vectores 

xhatsum = 0;
xhat2 = zeros(n,n*2);
% xhat1=Mw_f Mi_f Mc_f Ms_f;

for i=1:2*n    
Mw_f=xhat1(1,i);
Mi_f=xhat1(2,i);
Mc_f=xhat1(3,i);
Ms_f=xhat1(4,i);
 T_f=xhat1(5,i);
    
%% ----------------------------Model--------------------
%***********------------ Brix------------*********************************
Bx_ml_f=100*(Ms_f+Mi_f)/(Ms_f+Mi_f+Mw_f);
Bx_mg_f=100*(Ms_f+Mi_f+Mc_f)/(Ms_f+Mi_f+Mc_f+Mw_f);
%**************---------Purity--------*********************************
P_ml_f=Ms_f/(Ms_f+Mi_f);
P_mg_f=(Ms_f+Mc_f)/(Ms_f+Mi_f+Mc_f);
%***********---Concentration --***********************************
C_sol_f=Ms_f/Mw_f;
%******************--Correlations-----------*****************************
%----Density--------------
dens_feed_f=(1000+p.Bx_feed*(200+p.Bx_feed)/54)*(1-0.036*(p.T_feed-20)/(160-p.T_feed));     %kg/m3 
%----Specific heat capacity ----------
Cp_feed_f=(4186.8-29.7*p.Bx_feed+4.61*p.Bx_feed*p.P_feed+0.075*p.Bx_feed*p.T_feed);         %J/(kg C)
%----Enthalpy-----
h_feed_f=Cp_feed_f*p.T_feed;                                                                 %J/kg
%********------------Mother Liquor--***********************                                                                   % C
%-----Density--------------
dens_ml_f=(1000+Bx_ml_f*(200+Bx_ml_f)/54)*(1-0.036*(T_f-20)/(160-T_f));
%----Specific heat capacity ------------------------------
Cp_ml_f=(4186.8-29.7*Bx_ml_f+4.61*Bx_ml_f*P_ml_f+0.075*Bx_ml_f*T_f);                        %J/(kg C)
%----Enthalpy----------------------------
h_ml_f=Cp_ml_f*T_f;                                                                         %J/kg
%************---MAGMA: mOTHER LIQUOR +CRISTALES------**********************                                                                    % C
%-----Density--------------
dens_mg_f=(1000+Bx_mg_f*(200+Bx_mg_f)/54)*(1-0.036*(T_f-20)/(160-T_f));
%----Specific heat capacity -----------------
Cp_mg=4187-29.309*Bx_mg_f;                %J/(kg C)
%--Entalphy--------------------------------------------
h_mg_f=Cp_mg*T_f;                                  %J/kg
%*******-------------------Volume----------************************
MT_f=Mc_f+Ms_f+Mi_f+Mw_f;
 V_f=MT_f/dens_mg_f; % m3
%%**************-------Crystal------------------***************************
%----Specific heat capacity ----------
Cp_c_f=1163.2+3.488*T_f;           %kJ/(kg C)
%-----Entalpia---------------------------------
h_c_f=Cp_c_f*T_f;                    %J/kg
%--Specific latent heat-------------------
Lc=h_ml_f-h_c_f;                     %J/kg
%% -------------------Vapor Water enthalpy---------------------------
h_vap_f=2499980-24186*p.Pvac+(1891.1+106.1*p.Pvac)*T_f;
%% ---------------------Water enthalpy---------------------------------
h_w_f=2323.3+4106.7*p.Tw;
%% --------------------Crystales Conten------------------------------
cc_f=Mc_f/MT_f;
%% --------------------Heat Steam Flow------------------------------
Fvap_f=Mvapor_f/p.dens_cw;
Fcw_f=Fvap_f/p.alfa_vap;
Fhs_f=p.alfa_q*Fcw_f;
%%  ----------------------Heat -------------------------------------------
Lcw_f=2323.3+4106.7*p.Tw;   %  Entalpia del agua condensado
Qs_f=p.dens_cw*Fhs_f*Lcw_f; 
Lvap_f=h_vap_f-h_w_f;
%% -----------------------Non Linear Maping---------------------
%-------Water------- Kg
dMw_f=p.dens_feed*Flujo*(1-p.Bx_feed/100)-Mvapor_f;
%----Impurities----------Masa impurezas acumula en el tacho: Kg-----
dMi_f=p.dens_feed*Flujo*(p.Bx_feed/100)*(1-p.P_feed/100); 
 %---Crystals---------: Kg-----
dMc_f=cc_f*(p.dens_feed*Flujo-Mvapor_f)+p.alfa_c;
%----Sucrose: Kg--------------
dMs_f=p.dens_feed*Flujo*(p.Bx_feed/100)*(p.P_feed/100)-dMc_f; 
% ------Temperture---------------------------
dT_f=((p.W+Qs_f+p.dens_feed*Flujo*(h_feed_f-h_ml_f)-Mvapor_f*Lvap_f+dMc_f*Lc)/(MT_f*Cp_mg))/100; % C

 %-------------------------------------------------------------------------
F1=dMw_f;
F2=dMi_f;
F3=dMc_f;
F4=dMs_f;
F5=dT_f;


    system = [F1;F2;F3;F4;F5];
    xhat2(:,i) = xhat1(:,i) +system*ts;
    xhatsum = xhatsum + xhat2(:,i);
end

%% (3c) Se  combinan  los  vectores para  obtener  la  estimaci?n a priori del estado en el instante k
xhat2k = (1/(2*n))*xhatsum;
% % %%  (3d) Se  estima  la  covarianza a  priori, 
%   a?adiendoal  final  de  la ecuaci?n para tener en cuenta el ruido del proceso:
% Estimation of the a priori covariance of the estimation error
PP1sum = zeros(n,n);
for i=1:2*n
    PP1sum = PP1sum + (xhat2(:,i) - xhat2k)*(xhat2(:,i) - xhat2k)';
end
PP1 = (1/(2*n))*PP1sum + Qe;

%% (4a) Se eligen sigma-puntosasociados a las estimaciones a priori
%Sigma points generation: 2nd set

% sqPP = sqrt(n*PP)
sqPP2 =  (chol(n*PP1))';

for i=1:2*n
    if i<=n
        xhat3(:,i) = xhat2k + sqPP2(:,i);
    else
        xhat3(:,i) = xhat2k - sqPP2(:,i-n);
    end
end

%%  (4b) Se utiliza la ecuaci?n de medida no lineal para transformar los sigma-puntos vectores
yhatsum = 0;

Mw_xh3=xhat3(1,:);
Mi_xh3=xhat3(2,:);
Mc_xh3=xhat3(3,:);
Ms_xh3=xhat3(4,:);
T_xh3=xhat3(5,:);

for i=1:2*n
    
    Bx_mg_xh3=100*(Ms_xh3(1,i)+Mi_xh3(1,i)+Mc_xh3(1,i))/(Ms_xh3(1,i)+Mi_xh3(1,i)+Mc_xh3(1,i)+Mw_xh3(1,i));
    Bx_ml_xh3=100*(Ms_xh3(1,i)+Mi_xh3(1,i))/(Ms_xh3(1,i)+Mi_xh3(1,i)+Mw_xh3(1,i));
       MT_xh3=Mc_xh3(1,i)+Ms_xh3(1,i)+Mi_xh3(1,i)+Mw_xh3(1,i);
        V_xh3=MT_xh3/dens_mg_f; % m3
      
    yhat(1,i) =Bx_ml_xh3;
    yhat(2,i) =V_xh3;
    yhat(3,i) =T_xh3(i);  % Unable to perform assignment because the size of the left side is 1-by-1 and the size of the right side is 1-by-10.
    
    yhatsum = yhatsum + yhat(:,i);
    
%  yhat(i) = [0 0 0 0 1]*xhat3(:,i);
%  yhatsum = yhatsum + yhat(i);
end

%% (4c) Se   combinan   los   vectores para   obtener   las medidas predichas en el tiempo .
yhatk = (1/(2*n))*yhatsum;

%% (4d) Se  estima  la  covarianza  del  error,  a?adiendo para  tener  en cuenta el ruido de la medida:
% Covariance of the predicted measurement

Pysum = 0;
for i=1:2*n
%     Pysum = Pysum + (yhat(i) - yhatk)*(yhat(i) - yhatk)';
     Pysum = Pysum + (yhat(:,i) - yhatk)*(yhat(:,i) - yhatk)';
end

Py = (1/(2*n))*Pysum + Re;

%% (4e) Se estima la covarianza cruzada 
% Cross covariance between \hat{x}_{k}^{-} and \hat{y}_{k}

PPxysum = 0;

for i=1:2*n
    PPxysum = PPxysum + (xhat3(:,i) - xhat2k)*(yhat(:,i) - yhatk)';
end

Pxy = (1/(2*n))*PPxysum;

%% (4f) Actualization:
Kk = Pxy*pinv(Py);
xhatf = xhat2k + Kk*(y - yhatk);
PP = PP1 - Kk*Py*Kk';

xe = xhatf;
Pxe=PP;
ze = [xe; Pxe(:,1); Pxe(:,2); Pxe(:,3); Pxe(:,4); Pxe(:,5) ];
end
