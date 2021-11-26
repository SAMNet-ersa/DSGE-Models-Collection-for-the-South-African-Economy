%----------------------------------------------------------------
% The optimal monetary and macroprudential policies in an estimated DSGE Model for South Africa
%----------------------------------------------------------------

%----------------------------------------------------------------
% Defining Variables and Parameters
%----------------------------------------------------------------
var c c_s c_b c_e c_f u_cs u_cb u_ce u_cf h_s h_b h_e n_s n_b w_s w_b q d l_b l_e l r r_b r_e sprd_b sprd_e 
    zeta_b zeta_e y pi x kappa lambda_b lambda_e lambda_f a_j z gamma_b gamma_e 
    Lb_obs Le_obs NPLb_obs NPLe_obs pi_obs q_obs R_obs Y_obs; //varepsilon_b varepsilon_e

varexo ups_varepsilonb         //Household Borrowers repayment shock
       ups_varepsilone         //Entrepreneur repayment shock
       ups_j                   //Housing Demand shock
       ups_gammab              //Household LTV shock
       ups_gammae              //Entrepreneur LTV shock
       ups_z                   //Technology shock
       ups_r                   //Monetary policy shock
       ups_pi                  //Cost-push shock
       
       ups_zetab_me            //Entrepreneur repayment shock
       ups_zetae_me;           //Entrepreneur repayment shock
    
parameters beta_s beta_b beta_e beta_f eta j tau nu sigma theta iota_pi m_b m_e w_kappab w_kappae kappa_ss chi_yk chi_lk phi_r phi_pi phi_y rho_j rho_z rho_gammab rho_gammae 
           rho_zetab rho_zetae zeta_bss zeta_ess chi_zetab chi_zetae pi_ss X_ss phi_ds phi_df phi_le phi_ef phi_lb phi_bf vartheta_b vartheta_e
           rho_varepsilonb rho_varepsilone; 
            
%----------------------------------------------------------------
% Assigning parameters values based on the literature - Calibration (2)
% Without portfolio adj. costs in HHs and Entrep budget constraints
%----------------------------------------------------------------

beta_s = 0.995;
beta_b = 0.97;    //(0.94 - 0.98): See Iacoviello (2005), Gerali et al (2010), Angelini et al (2014), Gupta and Sun (2018), Rubio and Carrasco-Gallego (2016).
beta_e = 0.96;
beta_f = 0.955;   //In Rubio and Carrasco-Gallego (2016), this value is 0.965
eta = 0.50;       
j = 0.12;         //0.12 also works with m_b=0.8
nu = 0.070;       //Literature suggest values in this ragnge (0.05 - 0.12)
m_b = 0.80;  
m_e = 0.60; 
tau = 2;
sigma = 0.31;     //(0.28 - 0.34)
kappa_ss = 0.13;  //Historical average is 0.13 while the requirement is 0.105 
w_kappab = 1;
w_kappae = 1;

rho_j = 0.90;    //0.99 generates -ve effects on H-Savers and +ve effects on H-Borrowers (and correspond nicely with phi_if = 0.05), the opposite is true for 0.8 (and correspond nicely with phi_if = 0.5)
rho_z = 0.90;
rho_gammab = 0.7;
rho_gammae = 0.7;
rho_varepsilonb = 0.75;
rho_varepsilone = 0.75; 
rho_zetab = 0.70;
rho_zetae = 0.70;

zeta_bss = 0.04;   //Steady state NPLs, in data 0.35 but Steinbach calibrate it at 0.6 for NPL shock to have significant effect - check this with SARB
zeta_ess = 0.035;
chi_zetab = 0.75;   
chi_zetae = 0.75;   
vartheta_b  = 0.50;
vartheta_e  = 0.50;

theta = 0.75;
iota_pi = 0.5;
phi_r = 0.73;
phi_pi = 1.5;
phi_y  = 0.5; 
chi_lk = 0.00;  
chi_yk = 0.00;

phi_ds = 0.00;
phi_lb = 0.00;
phi_le = 0.00;
phi_df = 0.00;
phi_bf = 0.00;  
phi_ef = 0.00;  

pi_ss = 1.016; 
X_ss = 1.1;

%----------------------------------------------------------------
% Defining the model: first are composite parameters to handle parameter dependence correctly using (#-operator). 
%----------------------------------------------------------------

model(linear);
%-----------------------------
//Composite parameters
%-----------------------------
# R_ss = pi_ss/beta_s;
# lambda_fss = (beta_s - beta_f)/beta_s;
# R_bss = ((1 - lambda_fss*(1 - w_kappab*kappa_ss))/(beta_f*(1-zeta_bss) - lambda_fss*zeta_bss*(1-w_kappab*kappa_ss)))*pi_ss;
# R_ess = ((1 - lambda_fss*(1 - w_kappae*kappa_ss))/(beta_f*(1-zeta_ess) - lambda_fss*zeta_ess*(1-w_kappae*kappa_ss)))*pi_ss;
# sprd_bss = 4*(R_bss - R_ss);
# sprd_ess = 4*(R_ess - R_ss);

# F_b1 = beta_b*(1 - zeta_bss*(1 - vartheta_b))*(R_bss/pi_ss);
# F_b2 = beta_b + m_b*((pi_ss/R_bss) - beta_b*(1 - zeta_bss*(1 - vartheta_b)));

# F_e1 = beta_e*(1 - zeta_ess*(1 - vartheta_e))*(R_ess/pi_ss);
# F_e2 = beta_e + m_e*((pi_ss/R_ess) - beta_e*(1 - zeta_ess*(1 - vartheta_e)));

# F_fb1 = beta_f*(1 - zeta_bss)*(R_bss/pi_ss);
# F_fb2 = F_fb1 - (1 - w_kappab*kappa_ss)*lambda_fss*(R_bss/pi_ss)*zeta_bss;
# F_fb3 = (1 - zeta_bss*(R_bss/pi_ss));
# F_fe1 = beta_f*(1 - zeta_ess)*(R_ess/pi_ss);
# F_fe2 = F_fe1 - (1 - w_kappae*kappa_ss)*lambda_fss*(R_ess/pi_ss)*zeta_ess;
# F_fe3 = (1 - zeta_ess*(R_ess/pi_ss));

# z_1 = j/(1-beta_s);
# z_2 = j/(1-beta_b-m_b*((pi_ss/R_bss)-beta_b*(1 - zeta_bss*(1 - vartheta_b))));
# z_3 = 1/(1 - z_2*m_b*((pi_ss/R_bss) - (1 - zeta_bss*(1 - vartheta_b))));
# z_4 = (1 - beta_e - m_e*((pi_ss/R_ess)-beta_e*(1 - zeta_ess*(1 - vartheta_e))))/beta_e;
# z_5 = ((1-w_kappae*kappa_ss)*((R_ss/pi_ss) - 1)*(1 - zeta_ess*(R_ess/pi_ss)))*m_e*(pi_ss/R_ess)*(nu/z_4)*(1/X_ss) + (X_ss - 1)/X_ss;
# z_6 = (1-w_kappab*kappa_ss)*((R_ss/pi_ss) - 1)*(1 - zeta_bss*(R_bss/pi_ss))*(pi_ss/R_bss)*z_2*z_3*m_b;
# z_7 = ((1-sigma)*(1-nu))/X_ss;
# z_8 = 1 + (z_5/z_7) + z_6*(sigma/(1-sigma));
# z_9 = nu + ((pi_ss/R_ess) - (1 - zeta_ess*(1 - vartheta_e)))*(nu/z_4)*m_e;

# N_sss = 1/(1 + tau*z_8);
# N_bss = 1/(1 + tau*z_3);
# H_ess = nu/(nu + z_1*z_4*z_8*(1-sigma)*(1-nu) + z_2*z_3*z_4*sigma*(1-nu));
# H_sss = (z_1*z_4*z_8*(1-sigma)*(1-nu))/(nu + z_1*z_4*z_8*(1-sigma)*(1-nu) + z_2*z_3*z_4*sigma*(1-nu));
# H_bss = 1 - H_sss - H_ess;
# Y_ss = (H_ess^nu)*(((N_bss^sigma)*(N_sss^(1-sigma)))^(1-nu));
# qHs_Y = z_1*z_7*z_8;
# qHb_Y = z_2*z_3*sigma*(1-nu)*(1/X_ss);
# qHe_Y = (nu/z_4)*(1/X_ss);
# WsNs_Y = z_7;
# WbNb_Y =  sigma*(1-nu)*(1/X_ss);

# Lb_Y = z_2*z_3*m_b*sigma*(1-nu)*(pi_ss/R_bss)*(1/X_ss);
# L_bss = Lb_Y*Y_ss;
# Le_Y = m_e*(pi_ss/R_ess)*(nu/z_4)*(1/X_ss);
# L_ess = Le_Y*Y_ss;
# L_ss = L_bss + L_ess;

# D_Y = (1 - w_kappab*kappa_ss)*(1 - zeta_bss*(R_bss/pi_ss))*Lb_Y + (1 - w_kappae*kappa_ss)*(1 - zeta_ess*(R_ess/pi_ss))*Le_Y;

# Cs_Y = z_7*z_8;
# Cb_Y = z_3*sigma*(1-nu)*(1/X_ss);
# Ce_Y = z_9*(1/X_ss);
# Cf_Y = (1 - (R_ss/pi_ss))*D_Y - (1 - (R_bss/pi_ss)*(1 - zeta_bss))*Lb_Y - (1 - (R_ess/pi_ss)*(1 - zeta_ess))*Le_Y;
# C_Y = Cs_Y + Cb_Y + Ce_Y + Cf_Y;

%-----------------------------
//Model equations
%-----------------------------
//Patient Households
u_cs = - (1/(1 - eta))*(c_s - eta*c_s(-1));                                //Marginal ulity of consumption

u_cs(+1) - u_cs + r - pi(+1) = phi_ds*(d - d(-1));                         //FOC Deposits 

q = (1-beta_s)*(a_j - h_s) + beta_s*(u_cs(+1) + q(+1)) - u_cs;             //FOC Real Estate 

w_s = (N_sss/(1 - N_sss))*n_s - u_cs;                                      //FOC Labour 

Cs_Y*(c_s) = WsNs_Y*(w_s + n_s) - qHs_Y*(h_s - h_s(-1)) + ((X_ss-1)/X_ss)*y +(1/X_ss)*x - D_Y*(d - (R_ss/pi_ss)*(r(-1) - pi + d(-1)));   //Budget Constraint 

//Impatient Households
u_cb = - (1/(1 - eta))*(c_b - eta*c_b(-1));                                //Marginal ulity of consumption
 
l_b = q(+1) + h_b - (r_b - pi(+1)) + gamma_b;                              //Borrowing Constraint

F_b1*(u_cb(+1) - u_cb + r_b(+1) - pi(+1)) = (beta_b*zeta_bss*(1 - vartheta_b)*(R_bss/pi_ss))*zeta_b(+1) - phi_lb*(l_b - l_b(-1)) - (1 - F_b1)*(lambda_b - u_cb);  //FOC Loans
 
q = (1 - F_b2)*(a_j - h_b - u_cb) + F_b2*q(+1) + (F_b2 - beta_b)*(lambda_b - u_cb + gamma_b - r_b + pi(+1)) + beta_b*(u_cb(+1) - u_cb);    //FOC Real Estate

w_b = (N_bss/(1 - N_bss))*n_b - u_cb;                                      //FOC Labour 

Cb_Y*c_b = WbNb_Y*(w_b + n_b) - qHb_Y*(h_b - h_b(-1)) + Lb_Y*(l_b - ((R_bss/pi_ss)*(1-zeta_bss*(1 - vartheta_b)))*(l_b(-1) + r_b(-1) - pi)) + zeta_bss*(1 - vartheta_b)*(R_bss/pi_ss)*Lb_Y*zeta_b;  //Budget Constraint

zeta_b = rho_zetab*zeta_b(-1) - chi_zetab*(y - y(-1)) + ups_varepsilonb; 

//zeta_b = - chi_zetab*(y - y(-1)) + varepsilon_b;  //Also works

//zeta_b = rho_zetab*zeta_b(-1) - chi_zetab*y + ups_varepsilonb;

//Entrepreneurs
u_ce = - (1/(1 - eta))*(c_e - eta*c_e(-1));                                //Marginal ulity of consumption

y = z + nu*h_e(-1) + (1-nu)*sigma*n_b + (1-nu)*(1-sigma)*n_s;              //Entrep Production Function

l_e = q(+1) + h_e - (r_e(+1) - pi(+1)) + gamma_e;                          //Entrep Collateral Constraint 

F_e1*(u_ce(+1) - u_ce + r_e(+1) - pi(+1)) = beta_e*zeta_ess*(1 - vartheta_e)*(R_ess/pi_ss)*zeta_e(+1) - phi_le*(l_e - l_e(-1)) - (1 - F_e1)*(lambda_e  - u_ce); //FOC Loans

q = (1 - F_e2)*(y(+1) - x(+1) - h_e) + (1 + beta_e - F_e2)*(u_ce(+1) - u_ce) + F_e2*q(+1) + (F_e2 - beta_e)*(lambda_e - u_ce - r_e(+1) + pi(+1) + gamma_e);     //Entrep FOC Real Estate 

w_s = y - x - n_s;                                                         //Entrep FOC Labour_s 

w_b = y - x - n_b;                                                         //Entrep FOC Labour_b  

Ce_Y*c_e = Le_Y*(l_e - (R_ess/pi_ss)*(1-zeta_ess*(1 - vartheta_e))*(l_e(-1) + r_e - pi)) + zeta_ess*(1 - vartheta_e)*(R_ess/pi_ss)*Le_Y*zeta_e + (1/X_ss)*(y - x) - qHe_Y*(h_e - h_e(-1)) - WsNs_Y*(w_s + n_s) - WbNb_Y*(w_b + n_b);  //Entrep Budget Constraint  

zeta_e = rho_zetae*zeta_e(-1) - chi_zetae*(y - y(-1)) + ups_varepsilone; 

//zeta_e = - chi_zetae*(y - y(-1)) + varepsilon_e; //Also works

//zeta_e = rho_zetae*zeta_e(-1) - chi_zetae*y + ups_varepsilone;

//Bankers
u_cf = - (1/(1 - eta))*(c_f - eta*c_f(-1));                                //Marginal ulity of consumption

(beta_f*(R_ss/pi_ss))*(r - pi(+1)) = - (beta_f*(R_ss/pi_ss))*(u_cf(+1) - u_cf) - lambda_fss*(lambda_f - u_cf) - phi_df*(d - d(-1));  //Bank FOC Deposits 

F_fb2*(r_b - pi(+1)) = (beta_f*(R_bss/pi_ss) - F_fb2)*zeta_b(+1) - (1 - w_kappab*kappa_ss)*lambda_fss*F_fb3*(lambda_f - u_cf) + 
(w_kappab*kappa_ss*lambda_fss*F_fb3)*kappa - F_fb1*(u_cf(+1) - u_cf) + phi_bf*(l_b - l_b(-1));                                       //Bank FOC HH Lending

F_fe2*(r_e(+1) - pi(+1)) = (beta_f*(R_ess/pi_ss) - F_fe2)*zeta_e(+1) - (1 - w_kappae*kappa_ss)*lambda_fss*F_fe3*(lambda_f - u_cf) + 
(w_kappae*kappa_ss*lambda_fss*F_fe3)*kappa - F_fe1*(u_cf(+1) - u_cf) + phi_ef*(l_e - l_e(-1));                                       //Bank FOC Entrep Lending

D_Y*d = (1-w_kappab*kappa_ss)*Lb_Y*(F_fb3*l_b - zeta_bss*(R_bss/pi_ss)*(r_b - pi(+1) + zeta_b(+1))) - w_kappab*kappa_ss*F_fb3*Lb_Y*kappa
+  (1-w_kappae*kappa_ss)*Le_Y*(F_fe3*l_e - zeta_ess*(R_ess/pi_ss)*(r_e(+1) - pi(+1) + zeta_e(+1))) - w_kappae*kappa_ss*F_fe3*Le_Y*kappa;           //Bankers Capital Requirement/Collateral Constraint 

//Cf_Y*c_f = D_Y*(d - (R_ss/pi_ss)*(d(-1) + r(-1) - pi)) - Lb_Y*(l_b - ((1-zeta_bss)*(R_bss/pi_ss))*(l_b(-1) + r_b(-1) - pi)) - zeta_bss*(R_bss/pi_ss)*Lb_Y*zeta_b
//- Le_Y*(l_e - ((1 - zeta_ess)*(R_ess/pi_ss))*(l_e(-1) + r_e - pi)) - zeta_ess*(R_ess/pi_ss)*Le_Y*zeta_e;                                        //Bank Budget Constraint, but not needed  

//Housing Market
H_sss*h_s + H_bss*h_b + H_ess*h_e = 0;                                     //Total Real Estate Normalized to 1 

//Credit Market
L_bss*l_b + L_ess*l_e = L_ss*l;                                            //Total Credit 

//Interest rate spreads
sprd_b = 4*((R_bss/sprd_bss)*r_b - (R_ss/sprd_bss)*r);                     //Interest rate spread 

sprd_e = 4*((R_ess/sprd_ess)*r_e - (R_ss/sprd_ess)*r);                     //Interest rate spread 

//Total Consumption
C_Y*c = Cs_Y*c_s + Cb_Y*c_b + Ce_Y*c_e + Cf_Y*c_f;

y = Cs_Y*c_s + Cb_Y*c_b + Ce_Y*c_e + Cf_Y*c_f;
  
//Macroprudential Policies - Capital Requirements Rule
kappa = chi_lk*(l - y);

//Monetary Policy - Interest rate rule
r = phi_r*r(-1) + (1-phi_r)*((phi_pi*pi + phi_y*(y - y(-1)))) + ups_r;

//r = phi_r*r(-1) + (1-phi_r)*((phi_pi*pi(+1) + phi_y*(y - y(-1)))) + ups_r;

//Inflation dynamics 
pi = (iota_pi/(1 + beta_s*iota_pi))*pi(-1) + (beta_s/(1 + beta_s*iota_pi))*pi(+1) - (((1-theta)*(1-beta_s*theta))/(theta*(1 + beta_s*iota_pi)))*x + ups_pi;

//pi = beta_s*pi(+1) - (((1-theta)*(1-beta_s*theta))/theta)*x + ups_pi;                             

//Shock Processes
gamma_b = rho_gammab*gamma_b(-1) + ups_gammab;                             //Household LTV shock

gamma_e = rho_gammae*gamma_e(-1) + ups_gammae;                             //Entrep LTV shock

a_j = rho_j*a_j(-1) + ups_j;                                               //Housing Demand Shock

z = rho_z*z(-1) + ups_z;                                                   //Technology Shock

//varepsilon_b = rho_varepsilonb*varepsilon_b(-1) + ups_varepsilonb;         //Household borrower repayments shock (i.i.d shock)

//varepsilon_e = rho_varepsilone*varepsilon_e(-1) + ups_varepsilone;         //Entrepreneur repayment shock (i.i.d shock)

%----------------------------------------------------------------
% Measurement Equation
%---------------------------------------------------------------
Lb_obs = l_b - l_b(-1);
Le_obs = l_e - l_e(-1);
NPLb_obs = zeta_b + ups_zetab_me;
NPLe_obs = zeta_e + ups_zetae_me;
Y_obs = y - y(-1);
q_obs = q - q(-1);
pi_obs = pi;
R_obs = r;

end;

%---------------------------------------------------------------
%  Define shock variances
%---------------------------------------------------------------

shocks;
var ups_varepsilonb; stderr 0.15;
var ups_varepsilone; stderr 0.15;
var ups_j; stderr 0.03;
var ups_gammab; stderr 0.01;
var ups_gammae; stderr 0.01; 
var ups_z; stderr 0.01; 
var ups_r; stderr 0.01;
var ups_pi; stderr 0.01;

var ups_zetab_me; stderr 0.1;
var ups_zetae_me;  stderr 0.1;
end;

//resid;

steady;

check;

%---------------------------------------------------------------
% Model Estimation
%---------------------------------------------------------------

estimated_params;
//%                   Initial values,            PRIORS
eta                    , 0.50 , 0 ,  1  , beta_pdf     ,   0.50    ,  0.050 ;   
sigma                  , 0.29 , 0 ,  1  , beta_pdf     ,   0.30    ,  0.010 ;     
theta                  , 0.70 , 0 ,  1  , beta_pdf     ,   0.65    ,  0.020 ;  
iota_pi                , 0.50 , 0 ,  1  , beta_pdf     ,   0.50    ,  0.050 ;   
phi_r                  , 0.70 , 0 ,  1  , beta_pdf     ,   0.70    ,  0.050 ;  
phi_pi                 , 1.60 , 0 , Inf , gamma_pdf    ,   1.70    ,  0.050 ; 
phi_y                  , 0.50 , 0 , Inf , normal_pdf   ,   0.50    ,  0.050 ; 
//phi_ds                 , 0.05 , 0 , Inf , gamma_pdf    ,   0.05    ,  0.025 ;  
//phi_df                 , 0.05 , 0 , Inf , gamma_pdf    ,   0.05    ,  0.025 ;
//phi_lb                 , 0.05 , 0 , Inf , gamma_pdf    ,   0.05    ,  0.025 ; 
phi_bf                 , 0.05 , 0 , Inf , gamma_pdf    ,   0.25    ,  0.125;
//phi_le                 , 0.05 , 0 , Inf , gamma_pdf    ,   0.05    ,  0.025 ;
phi_ef                 , 0.05 , 0 , Inf , gamma_pdf    ,   0.25    ,  0.125 ;
chi_zetab              , 0.60 , 0 , Inf , gamma_pdf    ,   0.50    ,  0.100 ;
chi_zetae              , 0.60 , 0 , Inf , gamma_pdf    ,   0.50    ,  0.100 ;

vartheta_b             , 0.60 , 0 ,  1  , beta_pdf     ,   0.50    ,  0.100 ;
vartheta_e             , 0.60 , 0 ,  1  , beta_pdf     ,   0.50    ,  0.100 ;

rho_j                  , 0.80 , 0 ,  1  , beta_pdf     ,   0.70    ,  0.100 ;    
rho_z                  , 0.80 , 0 ,  1  , beta_pdf     ,   0.70    ,  0.100 ;     
rho_gammab             , 0.80 , 0 ,  1  , beta_pdf     ,   0.70    ,  0.100 ;     
rho_gammae             , 0.80 , 0 ,  1  , beta_pdf     ,   0.70    ,  0.100 ;
//rho_zetab              , 0.75 , 0 ,  1  , beta_pdf     ,   0.80    ,  0.100 ;     
//rho_zetae              , 0.75 , 0 ,  1  , beta_pdf     ,   0.80    ,  0.100 ;
//rho_varepsilonb        , 0.90 , 0 ,  1  , beta_pdf     ,   0.80    ,  0.100 ;     
//rho_varepsilone        , 0.90 , 0 ,  1  , beta_pdf     ,   0.80    ,  0.100 ;

stderr ups_varepsilonb , 0.10 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ; 
stderr ups_varepsilone , 0.10 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;
stderr ups_j           , 0.10 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;     
stderr ups_gammab      , 0.05 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;   
stderr ups_gammae      , 0.05 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;       
stderr ups_z           , 0.05 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ; 
stderr ups_r           , 0.05 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;    
stderr ups_pi          , 0.05 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;  
  
stderr ups_zetab_me    , 0.10 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;
stderr ups_zetae_me    , 0.10 , 0 , Inf , inv_gamma_pdf,   0.10   ,  0.250 ;

end;

%---------------------------------------------------------------
% Declare observables
%---------------------------------------------------------------

varobs Lb_obs Le_obs NPLb_obs NPLe_obs pi_obs q_obs R_obs Y_obs ;

identification(advanced=2);

//estimation(datafile=LM3_EstData_v8,smoother,bayesian_irf,irf=20,mh_conf_sig=0.95,mh_jscale=0.45,mode_compute=4,mh_drop=0.5,mode_check,prior_trunc=1e-32,mh_replic=250000,mh_nblocks=2,lik_init=2,graph_format=(eps,fig),moments_varendo,contemporaneous_correlation,conditional_variance_decomposition=[1,4,8,20,32],order=1) 

estimation(datafile=LM3_EstData_v8,endogenous_prior,smoother,bayesian_irf,irf=20,mh_conf_sig=0.95,mh_jscale=0.40,mode_compute=4,mh_drop=0.5,mode_check,prior_trunc=1e-32,mh_replic=250000,mh_nblocks=2,lik_init=1,graph_format=(eps,fig),moments_varendo,contemporaneous_correlation,conditional_variance_decomposition=[1,4,8,20],order=1) 
y c pi q l_b l_e r r_b r_e l d sprd_b sprd_e zeta_b zeta_e x c_s c_b c_e c_f h_s h_b h_e w_s w_b n_s n_b Lb_obs Le_obs NPLb_obs NPLe_obs pi_obs q_obs R_obs Y_obs;
//l d sprd_b sprd_e zeta_b zeta_e x c_s c_b c_e c_f h_s h_b h_e    

shock_groups(name=group1);
'LTV shocks' = ups_gammab, ups_gammae;
'NPL shocks' = ups_varepsilonb, ups_varepsilone;
'Monetary shock' = ups_r;
'Housing demand shock' = ups_j;
'Cost-push shock' = ups_pi;
'Technology shock' = ups_z;
'Measurement Error' = ups_zetab_me, ups_zetae_me;
end;

generate_trace_plots(2);

shock_decomposition(use_shock_groups=group1) Lb_obs Le_obs NPLb_obs NPLe_obs pi_obs q_obs R_obs Y_obs;

//plot_shock_decomposition(use_shock_groups=group1,write_xls) Lb_obs Le_obs NPLb_obs NPLe_obs pi_obs q_obs R_obs Y_obs;

