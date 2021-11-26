// Updated Main Estimation Code for HL2015 (Paper 2) 03 Nov 2015 --- JBF
// SAMPLE PERIOD 1984Q01 - 2014Q04 // log-difference of same dataset on extended series // Annualized rates vs Quart. rates
// NOTE: phi_w=0.8 and phi_k=0.8 prev. fixed for estimation.
// UPDATES (2021): phi_w & phi_k now estimated; phi_psi, theta_w and gamma_w estimated too (all parameters are identified).
//         Model specification change: FOC wage eqn replaced with HH BC because flex price equilibrium implies wage-setting eqn should be the HH FOC for labour.
// ESTIMATION: both #gamma_h's work for mode_compute=4; implied version used
//             This version excludes observed interbank rate: restrictive identification of bank capital req shocks and bank balance sheet adjustments on spreads       
//             measurement errors for obs data included, but negligible for all obs (< 0.000%). Priors also more ``agnostic'' and less restricted i.t.o std error

var b c_h i pi i_h q_psi lambda_h l_h w h nu_h mrs_w
    y x c_e i_e lambda_e v k_e nu_e l_e q_k   

    i_com k_B l omega_B k_BL S_com S_e S_h
    c  
    xi_z xi_i xi_b mu_e mu_h xi_p xi_t xi_psi

    y_obs pi_obs q_psi_obs b_obs l_h_obs l_e_obs i_obs i_com_obs i_h_obs i_e_obs; //

varexo epsilon_z epsilon_i epsilon_p epsilon_b epsilon_nu_h epsilon_nu_e epsilon_h epsilon_e epsilon_t epsilon_psi
       me_y me_p me_lh me_le me_b me_qpsi ; //

parameters beta_h phi R R_h zeta_psi Nu_h phi_w eta gamma
           R_e beta_e gamma_e kappa_v alpha Nu_e phi_k delta_e

           beta_R theta_R varepsilon_p gamma_p
           beta theta_w varepsilon_w gamma_w

           beta_B delta_B kappa_k tau Nu_B kappa_e varepsilon_ess kappa_h varepsilon_hss L_hL L_eL BL phi_psi

           kappa_i kappa_pi kappa_y 
           CY K_BY PsiY K_eY phi_B LY BY
 
           rho_z rho_i rho_b rho_e rho_h rho_nuh rho_nue rho_p rho_t rho_psi 
           pi_ss i_ss i_com_ss i_h_ss i_e_ss; // CY C_hY C_eY       

//Calibrated Parameters

//households
beta_h = 0.97;      // average of 0.99 and 0.95 is 0.97: see Iacoviello (2005), p.752
R_h = 1.02132;      // average of steady-state Rate of Returns from data (R_h + R)/2 = (1.03263427 + 1.01)/2 = 1.02132
//beta_h = 0.96;      //discount factor for borrower
//R_h = 1.0302;       // gross return to HH loans #Lambda = 1/R_h - beta_b must be bigger than zero in steady state. Therefore R_h corresponds to beta_b
phi = 0.65;          // habit formation
R = 1.01;           // 1.01187 gross real return to deposits
zeta_psi = 0.02988; // 1.02988^(0.25) = 1.0074 
Nu_h = 0.5;        // loan-to-value for HH
phi_w = 0.8;       // weight on wage income in borrowing constraint
eta = 1;            // inverse labour elasticity
gamma = 2;          // coefficient of relative risk aversion

//entrepreneurs
R_e = 1.03876;       // 1/1.03876 = 0.96269 gross return to Entrepreneurs loans see R_h above
beta_e = 0.954;     // entrepreneurs discount factor must be less than 1/R_e
gamma_e = 0.9;        // coefficient of relative risk aversion (see p.743, footnote 7, Iacoviello (2005))
kappa_v = 2;        // physical capital adjustment cost
delta_e = 0.025;    // depreciation rate of capital
alpha = 0.33;        // share of capital in firm production
Nu_e = 0.6;         // loan-to-value ratio for entrepreneurs
phi_k = 0.8;       // weight on the value of capital in borrowing constraint
varepsilon_p = 11;   // 9/8 is 12.5% markup: average of ---> prev. 7.667 - 21 (Iacoviello 2005) or 6 (Christiano et al 2010 for intermediate goods) price elasticity of demand(substitution) across differentiated retail goods

//retailers and unions

beta_R = 0.99;      // Retailers discount factor same as Saver household (0.99) or 1/R
theta_R = 0.8;     // probability of retail prices remaining unchanged in NKPC, theta_R = 0 for fully flexible prices
gamma_p = 0.5;      // degree of price indexation

beta = 0.97;        // Unions discount factor: equal to households
theta_w = 0.75;     // probability of wages remaining unchanged, theta_w = 0 for fully flexible wages
gamma_w = 0.3;      // degree of wage indexation, gamma_w = 0 no indexation and gamma_w = 1 full indexation
varepsilon_w = 5;   // wage elasticity of substitution across different types of HH

//banking sector

beta_B = 0.986;         // Banking discount factor is equal to the 1/(R_com) = 1/(1.0138375): a markup of 0.38375 over the 3MTB rate
delta_B = 0.1044;        // Gerali 0.1049 USmodel 0.1044 Cost for managing the banks capital position
//PsiK_B = 2.9853;       // 2.9853*0.3 compl with omegaK_B = 0.1044
tau = 0.11;              // target capital-to-loans ratio or capital requirement (leverage ratio) this is inconsistent with the data setup: should be 0.5-0.7
phi_psi = 0.25;           // conversion of equity to capital in cap.acc.eqn.

Nu_B = 0.4;                // value-at-risk constraint LTV for banks (see Woodford (2010) for discussion)
kappa_k = 2;               // parameter governing adjustment costs in banking
kappa_e = 5;               // parameter governing firm loan rate adjustment costs
kappa_h = 10;               // parameter governing HH loan rate adjustment costs
varepsilon_ess = 1.348;    // USdata 3.875766 steady state markup, varepsilon_ess/(varepsilon_ess - 1) is the markup on the loan rate to entrepreneurs
varepsilon_hss = 1.442;   // USdata 3.263427 steady state markup, varepsilon_hss/(varepsilon_hss - 1) is the markup on the loan rate to HHs
 
L_hL = 0.46;             // HH loans to total loans ratio
L_eL = 0.54;             // Entrepreneurial loans to total loans ratio
BL = 0.89;                // asset-loan ratio equal to 1 - tau or data 0.856
//Psi_BL = 0.275;         // 2.2*0.11 // 0.5206 0.5488  PsiY/LY*0.3 = 0.5488*0.3 = 0.1646 //(1.04) parameter in retained earnings equation equity-loan ratio
//Omega_BL = 0.0275;      // 0.34*0.11 // 0.02 0.04 loan:earnings ratio of 26 taken from entrep ratio 0.0375
phi_B = 0.3;              // bank equity:total equity ratio

//monetary policy and aggregates

kappa_i = 0.65;       // Taylor rule coefficient on i
kappa_pi = 1.5;       // Taylor rule coefficient on pi
kappa_y = 0.25;       // Taylor rule coefficient on y
//kappa_pich = 0.35;    //

CY = 0.655;          // data 0.6789 parameter ratios in aggregate resource eqn
//C_hY = 0.6364;        // hh consump-output ratio (see p.760 Iacoviello (2005))
//C_eY = 0.0426;         // entrep consump-output ratio (see p.760 Iacoviello (2005))

K_BY = 0.17;         // 0.165(K_BY=0.165)/(LY=1.5) gives 0.11 (K_B/L = tau)
LY = 1.55;             // 1.5 prev. 1 - data 1.5  
BY = 0.584;             // 0.6 (LY - K_BY) or see data for MA + deposits = 0.6
PsiY = 0.8;        // 0.7279 historical average 0.9 total equity to output ratio

K_eY = 10.7;            // 10.7 Capital-output ratio

//shock parameters

rho_z = 0.8;              // AR(1) parameter for Productivity shock
rho_i = 0.5;             // AR(1) parameter for MP shock
rho_b = 0.75;              // AR(1) parameter for Asset shock
rho_e = 0.75;              // AR(1) parameter for Entrep interest elasticity shock
rho_h = 0.75;              // AR(1) parameter for HH interest elasticity  shock

rho_nuh = 0.75;
rho_nue = 0.75;           // AR(1) parameter for LTV shock

rho_psi = 0.75;           // AR(1) parameter for Equity shock
rho_p = 0.5;             // AR(1) parameter for Cost-push shock
rho_t = 0.75;              // AR(1) parameter for Capital-asset shock

pi_ss = 0.00563654;       // inflation q-on-q steady-state; measurement eqn below converts to annualized
i_ss = 0.0385;            // fed funds i_ss. Quart 0.00942630329947646; Annualized: 0.0385
i_com_ss = 0.04234;      // steady state interest rate level Quart 0.010346878; Annualized: 0.04234
i_h_ss = 0.074975;        // steady state interest rate level Quart: 0.018189296; Ann. 0.074975
i_e_ss = 0.0810984;        // steady state interest rate level Quart. 0.019646241 ; Ann. 0.0810984

model(linear);

#aa = 1/(1 + ((1/R_h) - beta_h)*Nu_h*phi_w);     // 0.99, note saver ratio = 1   // in FOC borrower labour parameter, equal to (W/P)/MRS_b; (1/R_h - beta_b) must be positive = 0.01977, the combination of beta_b R_h and varepsilon_hss must be satisfied.
#gamma_h = 1/(1 - beta_h*R);                   // = 49.26 in FOC steady-state for deposits borrowers - utility share(a) of asset:consumption ratio
//#gamma_h = 0.856;                            // 0.856253 asset:consumption ratio from data
//#Gamma_psi = 1/(1 - ((1/R_h) - beta_h)*Nu_h*(1 - phi_w));     // Households
#Gamma_psi = ((1/R_h) - beta_h)*Nu_h*(1 - phi_w);       // with no equity in utility
//#gamma_psih = (1 - beta_h*(1 + zeta_psi) - ((1/R_h) - beta_h)*Nu_h*(1 - phi_w)); //  = 1/(1 - 0.9671 - 0.00636) = 37.68  // in FOC equity parameter on shock - utility share(1-a) of equity:consumption ratio
//#gamma_psih = (1/1.14);               // inverse of 1.14 used in Paper 1
//#Lambda_hss = (1/R_h - beta_h);       // = 0.00912 steady state FOC L_h                                

#Upsilon_k = ((1/R_e) - beta_e)*Nu_e*phi_k;     // 1 - 0.0076   Nu_e is LK_e in paper  // parameter in FOC K_e
//#Upsilon_u = 1/(1 - Upsilon_k); //
//#Lambda_ess = (1/R_e - beta_e);       // = 0.00811 steady state FOC L_e

#X = varepsilon_p/(varepsilon_p - 1);                        //steady state markup

#mu_w = varepsilon_w/(varepsilon_w - 1);                     // steady state wage mark-up
//#chi_h = (1/(mu_w*(1 + ((1/R_h) - beta_h)*Nu_h*phi_w)));     // parameter in wage setting equation (removed)
#phi_u = (1/(1 + beta));                 // wage-setting equation parameter
#phi_uu = (((1 - theta_w)*(1 - theta_w*beta))/((theta_w)*(1 + varepsilon_w*eta))); // wage-setting equation parameter

#HY = (1 - alpha)/X;
//#LK_e = Nu_e;     // note: assumed from paper, i.e. from borrowing constraint

#PsiK_B = (1 - delta_B); // 
#Psi_BL = PsiK_B*tau;
#Omega_BL = delta_B*tau;

#C_hY = (1 - alpha)/X - (R_h - 1)*L_hL*LY + (R - 1)*BY + zeta_psi*PsiY;           // now 0.6229 prev 0.65 hh consump-output ratio (see p.760 Iacoviello (2005))
//#C_eY = (alpha)/X - (R_e - 1)*L_eL*LY - delta_e*K_eY - zeta_psi*(1 - phi_B)*PsiY; // entrep consump-output ratio (see p.760 Iacoviello (2005))
//#CY = C_hY + C_eY;
#C_eY = CY - C_hY;  // approx 0.043

//Equations

//Households

//w(-1) = aa*(gamma/(1-phi))*(c_h - phi*c_h(-1)) + eta*h - (((1/R_h) - beta_h)*Nu_h*phi_w)*aa*(lambda_h + nu_h);      // FOC H with w(-1) see p. 20 dynare manual
((1/R_h) - beta_h)*lambda_h = beta_h*((gamma/(1-phi))*(c_h(+1) - phi*c_h) + pi(+1)) - (1/R_h)*((gamma/(1-phi))*(c_h - phi*c_h(-1)) + i_h);

b = (gamma_h)*(gamma/(1-phi))*(c_h - phi*c_h(-1)) + (beta_h*R*gamma_h)*(i - pi(+1) - (gamma/(1-phi))*(c_h(+1) - phi*c_h)) - xi_b;

q_psi = (q_psi(+1) - (gamma/(1-phi))*(c_h(+1) - phi*c_h)) + (gamma/(1-phi))*(1/(1 - Gamma_psi))*(c_h - phi*c_h(-1)) + (Gamma_psi/(1 - Gamma_psi))*(lambda_h + nu_h) - xi_psi;                                    // no equity in utility

//l_h = (phi_w/R_h)*(w(-1) + h) + ((1 - phi_w)/R_h)*q_psi - i_h + (1/R_h)*nu_h;        // borrowing constraint with w(-1) see p. 20 dynare manual
l_h = (phi_w/R_h)*(w + h) + ((1 - phi_w)/R_h)*q_psi - i_h + (1/R_h)*nu_h;   

C_hY*c_h = ((1 - alpha)/X)*(y - x) + BY*R*(i(-1) + b(-1) - pi) - BY*b 
        - L_hL*LY*R_h*(i_h(-1) + l_h(-1) - pi) + L_hL*LY*l_h
        + PsiY*(zeta_psi*q_psi - xi_psi); // HH flow of funds with h substituted out, p.761 Iacoviello (2005)

//Entrepreneurs

//h = y - x - w(-1);        // with w(-1) see p. 20 dynare manual
h = y - x - w;          // labour demand equation FOC H
((1/R_e) - beta_e)*lambda_e = beta_e*(gamma_e*(c_e(+1)) + pi(+1)) - (1/R_e)*(gamma_e*(c_e) + i_e); // FOC L_e

(1 - Upsilon_k)*(v - k_e(-1)) = (beta_e)*(v(+1) - k_e) + ((1 - beta_e*(1 - delta_e) - Upsilon_k)/(kappa_v))*(y(+1) - x(+1) - k_e) + ((Upsilon_k)/(kappa_v))*(lambda_e + nu_e) 
                            + ((beta_e*(1 - delta_e)*gamma_e)/(kappa_v))*(c_e - c_e(+1));  // FOC K_e determines the investment schedule
q_k = kappa_v*(v - k_e(-1)) - gamma_e*c_e;   // Shadow price of capital

y = alpha*k_e(-1) + (1 - alpha)*h + xi_z;        // production function
l_e = (phi_k/R_e)*(q_k + k_e(-1)) + ((1 - phi_k)/R_e)*(q_psi) - i_e + (1/R_e)*nu_e; // Borrowing constraint

k_e = (1 - delta_e)*k_e(-1) + delta_e*v;                 // Capital accumulation equation

C_eY*c_e = (alpha/X)*(y - x) - L_eL*LY*R_e*(i_e(-1) + l_e(-1) - pi) + L_eL*LY*l_e 
           - delta_e*K_eY*v - (1-phi_B)*PsiY*zeta_psi*(q_psi); //entrep flow of funds with h substituted out, p.761 Eq.(A9)Iacoviello (2005)

//Retailers
   // NKPC with indexation

pi = (beta_R/(1 + beta_R*gamma_p))*pi(+1) + (gamma_p/(1 + beta_R*gamma_p))*pi(-1) - (((1 - theta_R)*(1 - theta_R*beta_R))/((1 + beta_R*gamma_p)*theta_R))*x + xi_p;
   
//Unions and Wage setting equations
   // Wages with indexation

//w = (phi_u)*w(-1) + (phi_u*beta)*w(+1) + (phi_u*beta)*pi(+1) - (phi_u)*pi - (phi_u*theta_w*beta*gamma_w)*pi 
//     + (phi_u*gamma_w)*pi(-1) + (phi_uu)*((gamma/(1-phi))*(c_h - phi*c_h(-1)) + eta*h - w);

w = (phi_u)*w(-1) + (phi_u*beta)*w(+1) + (phi_u*beta)*pi(+1) - (phi_u)*pi - (phi_u*theta_w*beta*gamma_w)*pi 
     + (phi_u*gamma_w)*pi(-1) + (phi_uu)*(mrs_w - w);

mrs_w = aa*(gamma/(1-phi))*(c_h - phi*c_h(-1)) + eta*h - (((1/R_h) - beta_h)*Nu_h*phi_w)*aa*(lambda_h + nu_h);

//Banking Sector
 //Investment branch

k_B = (1 - delta_B)*k_B(-1) + delta_B*omega_B(-1) + phi_psi*(q_psi - q_psi(-1)) - (1 - phi_psi)*pi;   //  bank capital accumulation eqn

i_com = i - ((1/(R - 1))*kappa_k*((tau)^3))*(k_B - l - xi_t); // commercial bank loan rate with capital-asset shock    
S_com = i_com - i;                                            // spread between the commercial loan and policy rates
k_BL = k_B - l;                                               // bank capital-asset ratio

 //Commercial branch
i_e = ((kappa_e)/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*i_e(-1) + ((beta_B*kappa_e)/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*i_e(+1)
      + ((2*(varepsilon_ess - 1))/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*i_com 
      + (((1 - Nu_B)*(varepsilon_ess - 1))/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*mu_e;   // loan rate setting charged to entrepreneurs; the simplified case where varepsilon is a constant, drop last term. 

i_h = ((kappa_h)/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*i_h(-1) + ((beta_B*kappa_h)/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*i_h(+1)
      + ((2*(varepsilon_hss - 1))/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*i_com 
      + (((1 - Nu_B)*(varepsilon_hss - 1))/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*mu_h;   // loan rate setting charged to HHs; the simplified case where varepsilon is a constant, drop last term.

//i_e = (2/(1 - Nu_B))*i_com + mu_e;  //kappa_e = 0 for flexible adjustment
//i_h = (2/(1 - Nu_B))*i_com + mu_h;  //kappa_h = 0 for flexible adjustment

S_e = i_e - i_com;                      // spread between the entrepreneur loan and interbank rates
S_h = i_h - i_com;                      // spread between the household loan and interbank rates

 //Retained earnings net dividend payments
Omega_BL*(omega_B) = (R_h - 1)*L_hL*(i_h + l_h) + (R_e - 1)*L_eL*(i_e + l_e) - (R - 1)*BL*(i + b) - Psi_BL*zeta_psi*(q_psi);

//Monetary Policy
i = kappa_i*i(-1) + kappa_pi*(1 - kappa_i)*pi + kappa_y*(1 - kappa_i)*(y - y(-1)) + xi_i; // conventional nominal interest rate (Taylor-type) rule
//i = kappa_i*i(-1) + kappa_pi*(1 - kappa_i)*(pi(+1)) + kappa_pich*(1 - kappa_i)*(pi - pi(-1)) + kappa_y*(1 - kappa_i)*(y - y(-1)) + xi_i; // christiano2010-type nominal interest rate (Taylor-type) rule

//Aggregate resource constraints - Market clearing

y = C_hY*c_h + C_eY*c_e + delta_e*K_eY*v + delta_B*K_BY*k_B(-1);

CY*c = C_hY*c_h + C_eY*c_e;             // aggregate consumption
l = L_hL*l_h + L_eL*l_e;                // aggregate loans

//Shocks
xi_z = rho_z*xi_z(-1) + epsilon_z;                      // Productivity shock
xi_i = rho_i*xi_i(-1) + epsilon_i;                      // MP shock
xi_b = rho_b*xi_b(-1) + epsilon_b;                      // Asset shock

mu_e = rho_e*mu_e(-1) + epsilon_e;                      // Interest elasticity (markup) shocks
mu_h = rho_h*mu_h(-1) + epsilon_h;                      // Interest elasticity (markup) shocks

nu_h = rho_nuh*nu_h(-1) + epsilon_nu_h;                 // LTV shocks
nu_e = rho_nue*nu_e(-1) + epsilon_nu_e;                 // LTV shocks

xi_p = rho_p*xi_p(-1) + epsilon_p;                      // Cost-push shock
xi_t = rho_t*xi_t(-1) + epsilon_t;                      // capital-asset shock
xi_psi = rho_psi*xi_psi(-1) + epsilon_psi;              // Equity price shock

//Measurement equations

y_obs = y - y(-1) + 0.00389968 + me_y; //
pi_obs*4 = pi + pi_ss*4 + me_p;        // *4 annualized
q_psi_obs = q_psi - q_psi(-1) + 0.0144850805490961 + me_qpsi;  // 
b_obs = b - b(-1) + 0.00643461089818879 + me_b;  // 
l_h_obs = l_h - l_h(-1) + 0.00821458002965836 + me_lh;  // 
l_e_obs = l_e - l_e(-1) + 0.00618819398301963 + me_le;  // 

i_obs = i + i_ss;
i_com_obs = i_com + i_com_ss;
i_h_obs = i_h + i_h_ss;
i_e_obs = i_e + i_e_ss;

end;

//steady;
//check;

shocks;
%%% measurement errors
var me_y; stderr (0.006/(100/10));        // (0.026752018/(100/15)) sqrt
var me_p; stderr (0.0025/(100/10));        // (0.026752018/(100/15)) sqrt
var me_qpsi; stderr	(0.06/(100/15));     // (0.079364935/(100/25)) sqrt
var me_b; stderr (0.016/(100/15));        // 0.007 (0.026752018/(100/15)) sqrt
var me_lh; stderr (0.012/(100/10));        // 0.023613633 sqrt
var me_le; stderr (0.012/(100/10));        // 0.035875444 sqrt
end;

estimated_params;

  //Households
gamma, inv_gamma_pdf, 2, 0.5;          // coefficient of relative risk aversion
phi, beta_pdf, 0.65, 0.05;               // habit formation
//eta, gamma_pdf, 1, 0.5;
phi_w, beta_pdf, 0.8, 0.05;              // 0.45 weight in borrowing constraint
Nu_h, beta_pdf, 0.6, 0.05;               // works with 0.33 when phi_w is 0.05 0.65 LTV ratio

  // Retailers and Unions
theta_R, beta_pdf, 0.75, 0.05;           // price stickiness
gamma_p, beta_pdf, 0.5, 0.05;            // degree of wage indexation
//varepsilon_p, inv_gamma_pdf, 7.667, inf;
theta_w, beta_pdf, 0.75, 0.05;           // wage stickiness
gamma_w, beta_pdf, 0.5, 0.05;           // degree of wage indexation

   //Entrepreneurs
gamma_e, inv_gamma_pdf, 0.9, 0.1;        // coefficient of relative risk aversion
phi_k, beta_pdf, 0.8, 0.05;              // weight in borrowing constraint
Nu_e, beta_pdf, 0.6, 0.05;              // LTV ratio
phi_psi, gamma_pdf, 0.25, 0.05;          // equity-to-capital ratio

   //Banks
Nu_B, beta_pdf, 0.2, 0.05;             // value-at-risk LTV ratio

   //Adjustment costs
kappa_k, gamma_pdf, 2, 1;              // leverage deviation cost
kappa_e, gamma_pdf, 5, 2;              // Entrepreneur loan rate adjustment costs
kappa_h, gamma_pdf, 5, 2;              // HH loan rate adjustment costs
kappa_v, gamma_pdf, 2, 0.5;              // physical capital adjustment costs

  //Monetary Policy function
kappa_i, beta_pdf, 0.65, 0.05;          // Taylor rule coefficient on i
kappa_pi, gamma_pdf, 1.5, 0.05;        // Taylor rule coefficient on pi
kappa_y, beta_pdf, 0.25, 0.05;         // Taylor rule coefficient on y
//kappa_pich, beta_pdf, 0.35, 0.10;

   //Aggregate resource/flow of funds constraint
   
K_eY, gamma_pdf, 10.7, 0.2;              //
//CY, gamma_pdf, 0.6734, 0.1;             
//K_BY, gamma_pdf, 0.1804, 0.02;              
//LY, gamma_pdf, 0.9, 0.3;
//PsiY, gamma_pdf, 0.7672, 0.2; 
//BY, gamma_pdf, 0.8, 0.05;
//HY, gamma_pdf, 0.61, 0.05;           // i.e. (1-alpha)/X

//coefficients of AR(1) shocks

rho_z, beta_pdf, 0.75, 0.1;           // AR(1) parameter for productivity shock
rho_i, beta_pdf, 0.75, 0.1;           // AR(1) parameter for MP shock
rho_b, beta_pdf, 0.75, 0.1;            // AR(1) parameter for deposit shock
rho_e, beta_pdf, 0.75, 0.1;           // AR(1) parameter for entrep shock
rho_h, beta_pdf, 0.75, 0.1;           // AR(1) parameter for HH interest elasticity  shock

rho_nuh, beta_pdf, 0.75, 0.1;
rho_nue, beta_pdf, 0.75, 0.1;          // AR(1) parameter for LTV shock, added ad hoc

rho_psi, beta_pdf, 0.75, 0.1;          // AR(1) parameter for Equity shock
rho_p, beta_pdf, 0.75, 0.1;            // AR(1) parameter for Cost-push shock
rho_t, beta_pdf, 0.75, 0.1;           // AR(1) parameter for capital-asset shock shock

//innovations of the 10 AR(1) shocks

stderr epsilon_z, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_i, inv_gamma_pdf, 0.01, inf;
stderr epsilon_b, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_e, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_h, inv_gamma_pdf, 0.01, inf; 
  
stderr epsilon_nu_h, inv_gamma_pdf, 0.01, inf;
stderr epsilon_nu_e, inv_gamma_pdf, 0.01, inf;

stderr epsilon_psi, inv_gamma_pdf, 0.01, inf;
stderr epsilon_p, inv_gamma_pdf, 0.01, inf;
stderr epsilon_t, inv_gamma_pdf, 0.01, inf;

end;

// estimated_params_init(use_calibration);
// end;

varobs y_obs pi_obs q_psi_obs b_obs  l_h_obs l_e_obs  i_obs  i_h_obs i_e_obs; //    i_com_obs

// mode_file=Copy_of_HLbEst2015baseFULL_mode, mode_file=HLbEst2015baseFULL_mh_mode,(implied. gamma_h option w/ qpsi shock in HH BC) --> MAIN FILE

estimation(datafile=Estimation2015blndiff, mh_nblocks=3, mh_replic=100000, mh_drop=0.5, mh_jscale=0.28, // load_mh_file,
           mode_compute=0, mode_file=HLbEst2015baseFULL_mh_mode, mode_check, bayesian_irf, irf=20, moments_varendo, conditional_variance_decomposition=[1:20], 
           graph_format=fig, nodisplay) y pi i i_com i_h i_e S_com S_h S_e b l l_h l_e q_psi k_BL c c_h c_e k_e q_k k_B w h; // y_obs pi_obs q_psi_obs b_obs  l_h_obs l_e_obs;

// Shock Decomposition option

// shock_decomposition y_obs pi_obs S_com S_h S_e;

// //prior_mc=1 triggers the default point identification analysis; prior_mc>1 triggers the Monte Carlo mode;
// //identification(parameter set = prior mode);
//identification(parameter_set = posterior_mean, advanced=1, max_dim_cova_group=2);
