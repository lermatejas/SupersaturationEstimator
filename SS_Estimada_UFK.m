function SS_UFK = SS_Estimada_UFK(u)
  
Mw_est=u(1);
Mi_est=u(2);
Mc_est=u(3);
Ms_est=u(4);
T_est=u(5);

Bx_sat_est=64.447+0.08222*T_est+1.66169*10^(-3)*(T_est)^(2)-1.558*10^(-6)*(T_est)^(3)-4.63*10^(-8)*(T_est)^(4);
%----------------------------------------------
%a1=-0.063; c1=2.1 ;b1=0.935;                                             % Coeficientes de Lajos a1=-0.06265;b1=0.982;c1=2.1 para formula de Wiklund-Vavrinecz
 a1=-0.067; c1=2.8 ;b1=0.96;

q_NSW_est=Mi_est/Mw_est;                                                 % Relacion no azucares/agua ? Bx_ml/(100-Bx_ml) *(1-P_ml(t))
CSolb_mia_est=a1*q_NSW_est+b1+(1-b1)*exp(-c1*q_NSW_est);                 % Coeficiente Solubilidad, formula de Wiklund-Vavrinecz
SS_mia_est=(Ms_est/Mw_est)/((Bx_sat_est*CSolb_mia_est/(100-Bx_sat_est)));            % Sobresaturacion 
%--SS Teorica---------------
CSolb_ICUMSA_est=1-0.088*(Mi_est/Mw_est);
SS_ICUMSA_est=(Ms_est/Mw_est)/( Bx_sat_est/(100-Bx_sat_est))/CSolb_ICUMSA_est;
SS_UFK = SS_mia_est;
