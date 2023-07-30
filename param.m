function p=param()
%archivo de parametros
%%  Feed Flow 
p.Bx_feed=72;
p.P_feed=70;
p.T_feed=62;
p.dens_feed=(1000+p.Bx_feed*(200+p.Bx_feed)/54)*(1-0.036*(p.T_feed-20)/(160-p.T_feed));
p.Pol_feed=p.Bx_feed*p.P_feed;
%%  --CI FOOT ----------------------------------------------------
p.Bx_mg=88;                                          % B_mg    %
p.P_mg=75;  % p.P_mg=76.98                           % P_mg %
p.cc_o=26;p.cc_o=p.cc_o/100;                         % CC
p.Vpie=26.75;
p.Bx_ml=(100*p.cc_o-p.Bx_mg)/(p.cc_o-1);             % B_ml  
p.T0= 65; %                                          % Temperature
p.dens=(1000+p.Bx_mg*(200+p.Bx_mg)/54)*(1-0.036*(p.T0-20)/(160-p.T0)); % density kg/m3
%% --Condensate Flow  m3/seg-----------
p.Tw=62;                      % Temprature C     
p.dens_cw=1016.7-0.57*p.Tw;   % Density kg/m3 
%% ----Coef----
p.alfa_q=1.8177;
p.alfa_vap=1.6455;
p.alfa_c=0.9217;
%% Vacuum
p.Pvac=0.12;   % Bar
p.W=15;
end