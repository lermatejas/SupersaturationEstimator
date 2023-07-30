function  dxdt= PlantaBennet(Entrada)
p=param();

Mw=Entrada(1);   %Mass Water
Mi=Entrada(2);   %Mass Impurities 
Mc=Entrada(3);   %Mass Crystals
Ms=Entrada(4);   %Mass Sucrose
 T=Entrada(5);   %Temperature
Ff=Entrada(6);   %Feed Flow 
Mvap=Entrada(7); %Heat 

%% ---------------- Brix-------------------------------------------------
Bx_ml=100*(Ms+Mi)/(Ms+Mi+Mw);
Bx_mg=100*(Ms+Mi+Mc)/(Ms+Mi+Mc+Mw);
%% -----Purity ---------------
P_ml=Ms/(Ms+Mi);
P_mg=(Ms+Mc)/(Ms+Mi+Mc);
 %% -----Concentration---------------
C_sol=Ms/Mw;
 %% Supersaturation
Bx_sat=64.447+0.08222*T+1.66169*10^(-3)*(T)^(2)-1.558*10^(-6)*(T)^(3)-4.63*10^(-8)*(T)^(4);
%----------------------------------------------
%a1=-0.06265;b1=0.982;c1=2.1;                    % Lajos coef.        
a1=-0.067;  b1=0.96 ;c1=2.8 ;                    % Identified coefcent; Data form HM Sugar Mill
q_NSW=Mi/Mw;                                     % non_sugar/agua 
CSolb_mia=a1*q_NSW+b1+(1-b1)*exp(-c1*q_NSW);     % Solubility Eq. formula de Wiklund-Vavrinecz
SS_mia=(Ms/Mw)/((Bx_sat*CSolb_mia/(100-Bx_sat)));% Supersaturation 
%--SS Icumsa---------------
CSolb_ICUMSA=1-0.088*(Mi/Mw);
SS_ICUMSA=(Ms/Mw)/( Bx_sat/(100-Bx_sat))/CSolb_ICUMSA;
%% ---------------Feed--------------------------------------------
%----Density--------------
dens_feed=(1000+p.Bx_feed*(200+p.Bx_feed)/54)*(1-0.036*(p.T_feed-20)/(160-p.T_feed)); %kg/m3 
%----Specific enthalpy  ----------
Cp_feed=(4186.8-29.7*p.Bx_feed+4.61*p.Bx_feed*p.P_feed+0.075*p.Bx_feed*p.T_feed);     %J/(kg C)
%----enthalpy-----
h_feed=Cp_feed*p.T_feed;                                                              %J/kg
%% ---------------Mother Liquor-----------------------------
%-----Density--------------
dens_ml=(1000+Bx_ml*(200+Bx_ml)/54)*(1-0.036*(T-20)/(160-T));
%----Specific enthalpy ------------------------------
Cp_ml=(4186.8-29.7*Bx_ml+4.61*Bx_ml*P_ml+0.075*Bx_ml*T);                %J/(kg C)
%----enthalpy-----------------------------
h_ml=Cp_ml*T;                                                           %J/kg
%% ---------------MAGMA: mOTHER LIQUOR +CRISTALES--------------------------
%-----Density--------------
dens_mg=(1000+Bx_mg*(200+Bx_mg)/54)*(1-0.036*(T-20)/(160-T));
%----Specific enthalpy  -----------------
wc=Mc/(Mi+Ms+Mw+Mc); % fraccion de cristales
Cp_mg=4187-29.309*Bx_mg;                       %J/(kg C)
%-enthalpy --------------------------------------------
h_mg=Cp_mg*T;                                  %J/kg
%% Volume
MT=Mc+Ms+Mi+Mw;
V=MT/dens_mg; % m3
%% ------------------Crystal-----------------------------------------------
%----Specific enthalpy ----------
Cp_c=1163.2+3.488*T;           %kJ/(kg C)
%-----enthalpy---------------------------------
h_c=Cp_c*T;                    %J/kg
%--Specific latent heat-------------------
Lc=h_ml-h_c;                     %J/kg
%% -------------------Vapor Water enthalpy---------------------------
h_vap=2499980-24186*p.Pvac+(1891.1+106.1*p.Pvac)*T;
%% ---------------------Water enthalpy---------------------------------
h_w=2323.3+4106.7*p.Tw;
%% -------------------- Crystales Conten------------------------------
cc=Mc/MT;
%% --------------------Heat Steam Flow------------------------
%--Esto es par ameplear directamente 
Fvap=Mvap/p.dens_cw;
Fcw=Fvap/p.alfa_vap;
Fhs=p.alfa_q*Fcw;
%%  ----------------------Heat -------------------------------------------
Lcw=2323.3+4106.7*p.Tw;   
Qs=p.dens_cw*Fhs*Lcw; 
Lvap=h_vap-h_w;
%% -----------------------Mass Balance---------------------
%-------Water------- Kg
dMw=p.dens_feed*Ff*(1-p.Bx_feed/100)-Mvap;
 %---Crystals----------Masa cristales acumula en el tacho: Kg-----
dMc=cc*(p.dens_feed*Ff-Mvap)+p.alfa_c;
%----Sucrose Kg--------------
dMs=p.dens_feed*Ff*(p.Bx_feed/100)*(p.P_feed/100)-dMc; 
%----Impurities-----
dMi=p.dens_feed*Ff*(p.Bx_feed/100)*(1-p.P_feed/100); 
%% ---------------------------Energy---------------------------
dT=((p.W+Qs+p.dens_feed*Ff*(h_feed-h_ml)-Mvap*Lvap+dMc*Lc)/(MT*Cp_mg))/100;

%% Diferencials
dxdt=[
    dMw;
    dMi;
    dMc;
    dMs;
    dT;

    Bx_ml;
    Bx_mg;
    P_ml;
    P_mg;
    C_sol;
    SS_mia;
    SS_ICUMSA;
    MT;
    V;
    cc;
    ];
end