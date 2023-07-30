function SS_UFK = SS_Estimada_UFK(u)
  
Mw_est=u(1);
Mi_est=u(2);
Mc_est=u(3);
Ms_est=u(4);
T_real=u(5);
%--Saturation Brix 
Bx_sat_est=64.447+0.08222*T_real+1.66169*10^(-3)*(T_real)^(2)-1.558*10^(-6)*(T_real)^(3)-4.63*10^(-8)*(T_real)^(4);
%---- Coefficients------------------------------------------
 a1=-0.067; c1=2.8 ;b1=0.96; %Lajos coef. a1=-0.06265;b1=0.982;c1=2.1 
%-------------------------------
q_NSW_est=Mi_est/Mw_est;                                              % Non_sugar/agua
CSolb_est=a1*q_NSW_est+b1+(1-b1)*exp(-c1*q_NSW_est);                  % Solubility Coefficient, Wiklund-Vavrinecz Formula
SS_est=(Ms_est/Mw_est)/((Bx_sat_est*CSolb_est/(100-Bx_sat_est)));     % Estimated Supersaturation  
SS_UFK = SS_est;
