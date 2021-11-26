// Updated Main Estimation Code for HL2015 (Paper 2) 03 Nov 2015 --- JBF revise&resubmit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TAKEN OUT wage shock and replaced with capital-asset ratio shock xi_t
%%%%%%%%% Estimation Period: 1995Q1 - 2017Q2 (*South African data*)
%%%%%%%%% UPDATED estimatation file 12 obs. 13 shocks
%%%%%%%%% DATE STAMP: 2018/10/19
%%%%%%%%% THIS VERSION: k_BL_obs = xi_t + k_BL_ss; &  NO epsilon_k_BL bank capital shock (MAIN)
%%%%%%%%% THIS VERSION: includes output gap , based directly on (MAIN) w/ consistent priors & mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var b c_h i pi i_h q_psi lambda_h l_h w h nu_h
    y x c_e i_e lambda_e v k_e nu_e l_e q_k   
    mrs_w
    i_com k_B l omega_B k_BL S_com S_e S_h
    c  
    xi_z xi_w xi_v xi_g  xi_b mu_e mu_h xi_p xi_psi xi_t
    pi_wage winflationa inflationa interesta
    y_obs ygap_obs c_obs v_obs w_obs h_obs pi_obs q_psi_obs b_obs l_h_obs l_e_obs i_obs  k_BL_obs // ; //   xi_i  i_com_obs    i_h_obs i_e_obs

%%%%%% FLEX EQUILIBRIUM %%%%%%%
    ygap
    wn hn lambdan_h ln_h cn_h bn qn_psi //      
    yn vn kn_e cn_e qn_k lambdan_e ln_e // rn_k  xn 
    cn in in_com in_h in_e //          
    kn_B lsn omegan_B kn_BL Sn_com ; //  Sn_e Sn_h

%%%%%% FLEX EQUIL END %%%%%%%%
    
varexo epsilon_z epsilon_p epsilon_w epsilon_v epsilon_g epsilon_i epsilon_b  epsilon_t epsilon_nu_h epsilon_nu_e  epsilon_h epsilon_e epsilon_psi ; // epsilon_k_BL
       // me_y me_c me_v me_w me_h me_p me_lh me_le me_b me_qpsi ; //

parameters beta_h phi_h R R_h zeta_psi Nu_h phi_w eta_h gam
           R_e beta_e gamma_e kappa_v alph Nu_e phi_k delta_e

           beta_R theta_R varepsilon_p gamma_p
           beta theta_w varepsilon_w gamma_w

           beta_B delta_B kappa_k tau Nu_B kappa_e varepsilon_ess kappa_h varepsilon_hss L_hL L_eL BL phi_psi

            kappa_pi kappa_y 
            K_BY PsiY K_eY phi_B LY BY   // CY C_hY C_eY  
 
           rho_z rho_i rho_v rho_w rho_g rho_b rho_e rho_h rho_nuh rho_nue rho_p rho_t rho_psi cgy
           y_ss c_ss w_ss v_ss h_ss q_psi_ss b_ss l_h_ss l_e_ss pi_ss i_ss  i_h_ss i_e_ss k_BL_ss; //  i_com_ss  kappa_i

//Calibrated Parameters

//households
beta_h = 0.97;      // average of 0.99 and 0.95 is 0.97: see Iacoviello (2005), p.752
R_h = 1.018947;     // average of steady-state Rate of Returns from data (R_h + R)/2 = (1.018896 + 1.01)/2 = 
//beta_h = 0.96;      //discount factor for borrower
//R_h = 1.0302;       // gross return to HH loans #Lambda = 1/R_h - beta_b must be bigger than zero in steady state. Therefore R_h corresponds to beta_b
phi_h = 0.3;          // habit formation
R = 1.01;           // 1.01187 gross real return to deposits
zeta_psi = 0.019292; // 1.02988^(0.25) = 1.0074 // 0.008 used previously proportion of equity for dividends read stat that dividends are 15% of net earnings pre 2006, now down to 4%
Nu_h = 0.9;       // 0.8 loan-to-value for HH
phi_w = 0.85;       // 0.52 weight on wage income in borrowing constraint
eta_h = 2;            // inverse labour elasticity
gam = 1;          // coefficient of relative risk aversion

//entrepreneurs
R_e = 1.0131;       // 1/1.03876 = 0.96269 gross return to Entrepreneurs loans see R_h above
beta_e = 0.975;     // 0.975664 entrepreneurs discount factor must be less than 1/R_e
gamma_e = 1;        // 1.5 coefficient of relative risk aversion (see p.743, footnote 7, Iacoviello (2005))
kappa_v = 2;        // physical capital adjustment cost
delta_e = 0.025;    // depreciation rate of capital %%% HL=0.03 ... fix
alph = 0.3;        // 0.25 share of capital in firm production
Nu_e = 0.5;       // 0.8 loan-to-value ratio for entrepreneurs
phi_k = 0.3;       // 0. 5weight on the value of capital in borrowing constraint
varepsilon_p = 11;   // 9/8 is 12.5% markup: average of ---> prev. 7.667 - 21 (Iacoviello 2005) or 6 (Christiano et al 2010 for intermediate goods) price elasticity of demand(substitution) across differentiated retail goods

//retailers and unions

beta_R = 0.99;      // Retailers discount factor same as Saver household (0.99) or 1/R
theta_R = 0.85;     // probability of retail prices remaining unchanged in NKPC, theta_R = 0 for fully flexible prices
gamma_p = 0.25;      // degree of price indexation

beta = 0.97;        // Unions discount factor: equal to households
theta_w = 0.85;     // probability of wages remaining unchanged, theta_w = 0 for fully flexible wages
gamma_w = 0.5;      // degree of wage indexation, gamma_w = 0 no indexation and gamma_w = 1 full indexation
varepsilon_w = 5;   // wage elasticity of substitution across different types of HH

//banking sector

beta_B = 0.988;         // Banking discount factor is equal to the 1/(R_com) = 1/(1.012098): a markup of 0.38375 over the 3MTB rate
delta_B = 0.1044;        // Gerali 0.1049 USmodel 0.1044 Cost for managing the banks capital position
//PsiK_B = 2.9853;       // 2.9853*0.3 compl with omegaK_B = 0.1044
tau = 0.1;              // target capital-to-loans ratio or capital requirement (leverage ratio) this is inconsistent with the data setup: should be 0.5-0.7
phi_psi = 0.3;           // conversion of equity to capital in cap.acc.eqn.

Nu_B = 0.6;                // 0.25 value-at-risk constraint LTV for banks (see Woodford (2010) for discussion)
kappa_k = 2;               // 6 parameter governing adjustment costs in banking
kappa_e = 10;               // 1 parameter governing firm loan rate adjustment costs
kappa_h = 5;               // 2 parameter governing HH loan rate adjustment costs
varepsilon_ess = 8.206;    // 1.138767 quart // USdata 3.875766 steady state markup, varepsilon_ess/(varepsilon_ess - 1) is the markup on the loan rate to entrepreneurs
varepsilon_hss = 3.483;    // 1.402711 quart // USdata 3.263427 steady state markup, varepsilon_hss/(varepsilon_hss - 1) is the markup on the loan rate to HHs
 
L_hL = 0.52;             // HH loans to total loans ratio
L_eL = 0.48;             // Entrepreneurial loans to total loans ratio
BL = 0.89;                // asset-loan ratio equal to 1 - tau or data 0.856
//Psi_BL = 0.275;           // 2.2*0.11 // 0.5206 0.5488  PsiY/LY*0.3 = 0.5488*0.3 = 0.1646 //(1.04) parameter in retained earnings equation equity-loan ratio
//Omega_BL = 0.0275;        // 0.34*0.11 // 0.02 0.04 loan:earnings ratio of 26 taken from entrep ratio 0.0375
phi_B = 0.3;              // bank equity:total equity ratio

//monetary policy and aggregates

//kappa_i = 0.65;       // Taylor rule coefficient on i
kappa_pi = 1.5;       // Taylor rule coefficient on pi
kappa_y = 0.5;       // Taylor rule coefficient on y
//kappa_pich = 0.35;    //

//CY = 0.612;          // data 0.61 (model implied 0.7255) parameter ratios in aggregate resource eqn
//C_hY = 0.6364;        // hh consump-output ratio (see p.760 Iacoviello (2005))
//C_eY = 0.0426;        // entrep consump-output ratio (see p.760 Iacoviello (2005))

K_BY = 0.17;         // 0.165(K_BY=0.165)/(LY=1.5) gives 0.11 (K_B/L = tau)
LY = 0.66;             // 1.5 prev. 1 - data 1.5  
BY = 0.53;             // 0.6 (LY - K_BY) or see data for MA + deposits = 0.6
PsiY = 2;        //  historical average  total equity to output ratio

K_eY = 7.36;            // RSA V/Y = 0.184 // U.S. 10.7 Capital-output ratio

//shock parameters

rho_z = 0.9;              // AR(1) parameter for Productivity shock
rho_i = 0.85;             // AR(1) parameter for MP shock
rho_b = 0.85;              // AR(1) parameter for Asset shock
rho_e = 0.8;              // AR(1) parameter for Entrep interest elasticity shock
rho_h = 0.6;              // AR(1) parameter for HH interest elasticity  shock

rho_nuh = 0.7;           // AR(1) parameter for LTV shock
rho_nue = 0.75;           // AR(1) parameter for LTV shock

rho_psi = 0.75;           // AR(1) parameter for Equity shock
rho_p = 0.6;             // AR(1) parameter for Cost-push shock
rho_w = 0.25;             // AR(1) parameter for wage shock
rho_v = 0.8;              // AR(1) parameter for investment shock
rho_g = 0.815;           // AR(1) parameter for exogenous spending V shock
rho_t = 0.5;              // AR(1) parameter for Capital-asset shock

cgy = 0.3;               // tech shock corrrelation with exog spending

y_ss = 0.002393097;
pi_ss = 0.013976856;     // inflation steady-state
c_ss = 0.003727795;
v_ss = 0.006615576;
w_ss = -0.000858784;
h_ss = 0.001353938;
//g_ss = 0.002312112;

q_psi_ss = 0.008784239; // CHECK SERIES
b_ss = 0.006020798;     // CHECK SERIES
l_h_ss = 0.003415916;
l_e_ss = 0.0094405963;

i_ss = 0.022398;                  // policy i_ss       
//i_com_ss = 0.0245115952;        // OPTION BM: steady state Interbank NCD 12-month interest rate level
//i_com_ss = 0.02523945;            // OPTION TB: steady state 10-yr TB interest rate level
i_h_ss = 0.02748511;                // steady state HH (effective interest income) interest rate level
i_e_ss = 0.02523945;              // steady state Firm (10-yr TB) interest rate level
k_BL_ss = 0.0080117;             // steady-state Tier 1 bank capital-asset ratio

model(linear);

#aa = 1/(1 + ((1/R_h) - beta_h)*Nu_h*phi_w);     // 0.99, note saver ratio = 1   // in FOC borrower labour parameter, equal to (W/P)/MRS_b; (1/R_h - beta_b) must be positive = 0.01977, the combination of beta_b R_h and varepsilon_hss must be satisfied.
#gamma_h = 1/(1 - beta_h*R);                   // = 20.82 in FOC steady-state for deposits borrowers - utility share(a) of asset:consumption ratio
//#gamma_h = BY/CY;                                // 0.866 asset:consumption ratio from data
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
#phi_uu = (((1 - theta_w)*(1 - theta_w*beta))/((theta_w)*(1 + varepsilon_w*eta_h))); // wage-setting equation parameter

#HY = (1 - alph)/X;
//#LK_e = Nu_e;     // note: theoretically assumed from paper, i.e. from borrowing constraint

#PsiK_B = (1 - delta_B); // 
#Psi_BL = PsiK_B*tau;
#Omega_BL = delta_B*tau;

#C_hY = (1 - alph)/X - (R_h - 1)*L_hL*LY + (R - 1)*BY + zeta_psi*PsiY;     // SA: 0.61 (total cons) // 0.6737 now 0.6229 prev 0.65 hh consump-output ratio (see p.760 Iacoviello (2005))
#C_eY = (alph)/X - (R_e - 1)*L_eL*LY - delta_e*K_eY - zeta_psi*(1 - phi_B)*PsiY; // entrep consump-output ratio (see p.760 Iacoviello (2005))
#CY = C_hY + C_eY;
//#C_eY = CY - C_hY;  // approx 0.043

//Equations

%%%%%%%%%%% FLEX EQUILIBRIUM START %%%%%%%%%%%%% NEW (23 variable // flex price-flex (retail) rate) OLD (8 variables)
%%%% HH %%%%
wn = (gam/(1-phi_h))*(cn_h - phi_h*cn_h(-1)) + eta_h*hn - (((1/R_h) - beta_h)*Nu_h*phi_w)*aa*(lambdan_h + nu_h);  // *(lambdan_h + nu_h ); 
((1/R_h) - beta_h)*lambdan_h = beta_h*((gam/(1-phi_h))*(cn_h(+1) - phi_h*cn_h)) - (1/R_h)*((gam/(1-phi_h))*(cn_h - phi_h*cn_h(-1)) + in_h);
bn = (gamma_h)*(gam/(1-phi_h))*(cn_h - phi_h*cn_h(-1)) + (beta_h*R*gamma_h)*(in - (gam/(1-phi_h))*(cn_h(+1) - phi_h*cn_h)) - xi_b;  // 
ln_h = (phi_w/R_h)*(wn + hn) + ((1 - phi_w)/R_h)*qn_psi - in_h + (1/R_h)*nu_h;    // 
qn_psi = (qn_psi(+1) - (gam/(1-phi_h))*(cn_h(+1) - phi_h*cn_h)) + (gam/(1-phi_h))*(1/(1 - Gamma_psi))*(cn_h - phi_h*cn_h(-1)) + (Gamma_psi/(1 - Gamma_psi))*(lambdan_h + nu_h) - xi_psi;                                    // no equity in utility

CY*cn_h = ((1 - alph)/X)*(yn) + BY*R*(in(-1) + bn(-1)) - BY*bn 
        - L_hL*LY*R_h*(in_h(-1) + ln_h(-1)) + L_hL*LY*ln_h
        + PsiY*(zeta_psi*qn_psi - xi_psi); //HH flow of funds with h substituted out, p.761 Iacoviello (2005)

%%%% FIRMS %%%%
wn + hn = yn;         
 ((1/R_e) - beta_e)*lambdan_e = beta_e*(gamma_e*(cn_e(+1))) - (1/R_e)*(gamma_e*(cn_e) + in_e); 
 (1 - Upsilon_k)*(vn - kn_e(-1)) = (beta_e)*(vn(+1) - kn_e) + ((1 - beta_e*(1 - delta_e) - Upsilon_k)/(kappa_v))*(yn(+1) - kn_e) + ((Upsilon_k)/(kappa_v))*(lambdan_e + nu_e)   
                            + ((beta_e*(1 - delta_e)*gamma_e)/(kappa_v))*(cn_e - cn_e(+1)) + xi_v ; // 

// qn_k = kappa_v*(vn(+1) - kn_e) - gamma_e*cn_e(+1);
 qn_k = kappa_v*(vn - kn_e(-1)) - gamma_e*cn_e;
%% (1 - Upsilon_k)*(vn(-1) - kn_e(-1)) = (beta_e)*(vn - kn_e) + ((1 - beta_e*(1 - delta_e) - Upsilon_k)/(kappa_v))*(yn - kn_e) + ((Upsilon_k)/(kappa_v))*(lambdan_e(-1) + nu_e)   
%%                            + ((beta_e*(1 - delta_e)*gamma_e)/(kappa_v))*(cn_e(-1) - cn_e) + xi_v;  
%% qn_k = kappa_v*(vn - kn_e) - gamma_e*cn_e ;
// in_e = yn - kn_e(-1); // arbitrage condition on the real return to capital
// in_e(+1) = yn(+1) - kn_e;
// in = rho_i*in(-1) + kappa_y*(1 - rho_i)*(yn - yn(-1));

yn = alph*kn_e(-1) + (1 - alph)*hn + xi_z;  
ln_e = (phi_k/R_e)*(qn_k + kn_e(-1)) + ((1 - phi_k)/R_e)*(qn_psi) - in_e + (1/R_e)*nu_e ;    //   
//kn_e(+1) = (1 - delta_e)*kn_e + delta_e*vn;   
kn_e = (1 - delta_e)*kn_e(-1) + delta_e*vn; 
C_eY*cn_e = (alph/X)*(yn) - L_eL*LY*R_e*(in_e(-1) + ln_e(-1)) + L_eL*LY*ln_e 
           - delta_e*K_eY*vn - (1-phi_B)*PsiY*zeta_psi*(qn_psi); 

%%%% BANKS %%%%
//kn_B(+1) = (1 - delta_B)*kn_B + delta_B*omegan_B + phi_psi*(qn_psi(+1) - qn_psi);   
kn_B = (1 - delta_B)*kn_B(-1) + delta_B*omegan_B(-1) + phi_psi*(qn_psi - qn_psi(-1));   
in_com = in - ((1/(R - 1))*kappa_k*((tau)^3))*(kn_B - lsn - xi_t);  //     
Sn_com = in_com - in;                                            
kn_BL = kn_B - lsn; 
                                             
//in_e = in_com + mu_e ;  //                    
//in_h = in_com + mu_h ;  //   
// in_e = (2/(1 - Nu_B))*in_com + mu_e;  //kappa_e = 0 for flexible adjustment
// in_h = (2/(1 - Nu_B))*in_com + mu_h;  //kappa_h = 0 for flexible adjustment
in_e = ((kappa_e)/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*in_e(-1) + ((beta_B*kappa_e)/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*in_e(+1)
      + ((2*(varepsilon_ess - 1))/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*in_com 
      + (((1 - Nu_B)*(varepsilon_ess - 1))/((1 - Nu_B)*(varepsilon_ess - 1) + kappa_e + kappa_e*beta_B))*mu_e;   // loan rate setting charged to entrepreneurs; the simplified case where varepsilon is a constant, drop last term. 

in_h = ((kappa_h)/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*in_h(-1) + ((beta_B*kappa_h)/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*in_h(+1)
      + ((2*(varepsilon_hss - 1))/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*in_com 
      + (((1 - Nu_B)*(varepsilon_hss - 1))/((1 - Nu_B)*(varepsilon_hss - 1) + kappa_h + kappa_h*beta_B))*mu_h;   // loan rate setting charged to HHs; the simplified case where varepsilon is a constant, drop last term.
          
          
Omega_BL*(omegan_B) = (R_h - 1)*L_hL*(in_h + ln_h) + (R_e - 1)*L_eL*(in_e + ln_e) - (R - 1)*BL*(in + bn) - Psi_BL*zeta_psi*(qn_psi);

%%%% AGG %%%%
//yn = C_hY*cn_h + C_eY*cn_e + delta_e*K_eY*vn + xi_g;  // + delta_B*K_BY*kn_B
yn = CY*cn + delta_e*K_eY*vn + delta_B*K_BY*kn_B(-1) + xi_g;  // 
CY*cn = C_hY*cn_h + C_eY*cn_e;             
lsn = L_hL*ln_h + L_eL*ln_e; 

%%%%%%%%%%% FLEX EQUILIBRIUM END %%%%%%%%%%%%%

%%%%%%%%%%% FRICTIONS EQUILIBRIUM START %%%%%%%%%%%%%

//Households

//w = aa*(gam/(1-phi_h))*(c_h - phi_h*c_h(-1)) + eta_h*h - (((1/R_h) - beta_h)*Nu_h*phi_w)*aa*(lambda_h + nu_h);      // FOC H with w(-1) see p. 20 dynare manual
((1/R_h) - beta_h)*lambda_h = beta_h*((gam/(1-phi_h))*(c_h(+1) - phi_h*c_h) + pi(+1)) - (1/R_h)*((gam/(1-phi_h))*(c_h - phi_h*c_h(-1)) + i_h);

b = (gamma_h)*(gam/(1-phi_h))*(c_h - phi_h*c_h(-1)) + (beta_h*R*gamma_h)*(i - pi(+1) - (gam/(1-phi_h))*(c_h(+1) - phi_h*c_h)) - xi_b;

q_psi = (q_psi(+1) - (gam/(1-phi_h))*(c_h(+1) - phi_h*c_h)) + (gam/(1-phi_h))*(1/(1 - Gamma_psi))*(c_h - phi_h*c_h(-1)) + (Gamma_psi/(1 - Gamma_psi))*(lambda_h + nu_h) - xi_psi;                                    // no equity in utility

l_h = (phi_w/R_h)*(w + h) + ((1 - phi_w)/R_h)*q_psi - i_h + (1/R_h)*nu_h;        // borrowing constraint with w(-1) see p. 20 dynare manual

CY*c_h = ((1 - alph)/X)*(y - x) + BY*R*(i(-1) + b(-1) - pi) - BY*b 
        - L_hL*LY*R_h*(i_h(-1) + l_h(-1) - pi) + L_hL*LY*l_h
        + PsiY*(zeta_psi*q_psi - xi_psi); //HH flow of funds with h substituted out, p.761 Iacoviello (2005)

//Entrepreneurs

//h = y - x - w;          // labour demand equation FOC H
h = y - x - w;        // with w(-1) see p. 20 dynare manual
((1/R_e) - beta_e)*lambda_e = beta_e*(gamma_e*(c_e(+1)) + pi(+1)) - (1/R_e)*(gamma_e*(c_e) + i_e); // FOC L_e

%%%%%%%%
 (1 - Upsilon_k)*(v - k_e(-1)) = (beta_e)*(v(+1) - k_e) + ((1 - beta_e*(1 - delta_e) - Upsilon_k)/(kappa_v))*(y(+1) - x(+1) - k_e) + ((Upsilon_k)/(kappa_v))*(lambda_e + nu_e) 
                            + ((beta_e*(1 - delta_e)*gamma_e)/(kappa_v))*(c_e - c_e(+1)) + xi_v;  // FOC K_e determines the investment schedule
 q_k = kappa_v*(v - k_e(-1)) - gamma_e*c_e;   // Shadow price of capital
%%%%%%%% 
%%%%%%%% Investment Schedule timing of v kicked back as in GHW
%%
%% (1 - Upsilon_k)*(v(-1) - k_e(-1)) = (beta_e)*(v - k_e) + ((1 - beta_e*(1 - delta_e) - Upsilon_k)/(kappa_v))*(y - x - k_e) + ((Upsilon_k)/(kappa_v))*(lambda_e(-1) + nu_e) 
%%                            + ((beta_e*(1 - delta_e)*gamma_e)/(kappa_v))*(c_e(-1) - c_e) + xi_v ;  // FOC K_e determines the investment schedule
%% q_k = kappa_v*(v(-1) - k_e(-1)) - gamma_e*c_e(-1) ;   // Shadow price of capital   
%%%%%

y = alph*k_e(-1) + (1 - alph)*h + xi_z;        // production function
l_e = (phi_k/R_e)*(q_k + k_e(-1)) + ((1 - phi_k)/R_e)*(q_psi) - i_e + (1/R_e)*nu_e; // Borrowing constraint

//k_e(+1) = (1 - delta_e)*k_e + delta_e*v;                 // + delta_e*kappa_v*xi_v  // Capital accumulation equation
k_e = (1 - delta_e)*k_e(-1) + delta_e*v;                 // + delta_e*kappa_v*xi_v  // Capital accumulation equation

C_eY*c_e = (alph/X)*(y - x) - L_eL*LY*R_e*(i_e(-1) + l_e(-1) - pi) + L_eL*LY*l_e 
           - delta_e*K_eY*v - (1-phi_B)*PsiY*zeta_psi*(q_psi); //entrep flow of funds with h substituted out, p.761 Eq.(A9)Iacoviello (2005)

//Retailers
   // NKPC with indexation

pi = (beta_R/(1 + beta_R*gamma_p))*pi(+1) + (gamma_p/(1 + beta_R*gamma_p))*pi(-1) - (((1 - theta_R)*(1 - theta_R*beta_R))/((1 + beta_R*gamma_p)*theta_R))*x + xi_p;
   
//Unions and Wage setting equations // Wage inflation or Real wage setting
   // Wages with indexation

%%%%%%%%%% 
%% w = (phi_u)*w(-1) + (phi_u*beta)*w(+1) + (phi_u*beta)*pi(+1) - (phi_u)*pi - (phi_u*theta_w*beta*gamma_w)*pi 
%%     + (phi_u*gamma_w)*pi(-1) + (phi_uu)*((gam/(1-phi_h))*(c_h - phi_h*c_h(-1)) + eta_h*h - w) + xi_w;
%%%%%%%%%%

pi_wage - gamma_w*pi(-1) = beta*pi_wage(+1) - theta_w*beta*gamma_w*pi 
                           + phi_uu*(mrs_w - w) + xi_w; // + epsil_w 
%%%%%%%%%%
 // mrs_w = (gam/(1-phi_h))*(c_h - phi_h*c_h(-1)) + eta_h*h; // marginal rate of substitution in domestic wage inflation equation
 mrs_w = aa*(gam/(1-phi_h))*(c_h - phi_h*c_h(-1)) + eta_h*h - (((1/R_h) - beta_h)*Nu_h*phi_w)*aa*(lambda_h + nu_h);
%%%%%%%%%%

//Banking Sector
 //Investment branch

k_B = (1 - delta_B)*k_B(-1) + delta_B*omega_B(-1) + phi_psi*(q_psi - q_psi(-1)) - (1 - phi_psi)*pi ;   // - epsilon_k_BL bank capital accumulation eqn
//k_B(+1) = (1 - delta_B)*k_B + delta_B*omega_B + phi_psi*(q_psi(+1) - q_psi) - (1 - phi_psi)*pi(+1) ;   // - epsilon_k_BL bank capital accumulation eqn

i_com = i - ((1/(R - 1))*kappa_k*((tau)^3))*(k_B - l - xi_t);       // commercial bank loan rate with capital-asset shock    
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
i = rho_i*i(-1) + kappa_pi*(1 - rho_i)*pi + kappa_y*(1 - rho_i)*(y - y(-1)) + epsilon_i; // conventional nominal interest rate (Taylor-type) rule
//i = kappa_i*i(-1) + kappa_pi*(1 - kappa_i)*pi + kappa_y*(1 - kappa_i)*(y - y(-1)) + xi_i; // conventional nominal interest rate (Taylor-type) rule
//i = kappa_i*i(-1) + kappa_pi*(1 - kappa_i)*(pi(+1)) + kappa_pich*(1 - kappa_i)*(pi - pi(-1)) + kappa_y*(1 - kappa_i)*(y - y(-1)) + xi_i; // christiano2010-type nominal interest rate (Taylor-type) rule

//Aggregate resource constraints - Market clearing

y = C_hY*c_h + C_eY*c_e + delta_e*K_eY*v + delta_B*K_BY*k_B(-1) + xi_g;  // 

CY*c = C_hY*c_h + C_eY*c_e;             // aggregate consumption
l = L_hL*l_h + L_eL*l_e;                // aggregate loans

//Shocks
xi_z = rho_z*xi_z(-1) + epsilon_z;                      // Productivity shock
//xi_i = rho_i*xi_i(-1) + epsilon_i;                      // AR(1) MP shock // use i.i.d
xi_b = rho_b*xi_b(-1) + epsilon_b;                      // Asset shock

xi_w = rho_w*xi_w(-1) + epsilon_w;                      // (3) domestic wage markup shock
xi_v = rho_v*xi_v(-1) + epsilon_v;                      // (4) domestic investment shock
xi_g = rho_g*xi_g(-1) + epsilon_g + cgy*epsilon_z;      // (5) spending shock (a la SW (2007))

mu_e = rho_e*mu_e(-1) + epsilon_e;                      // Interest elasticity (markup) shocks
mu_h = rho_h*mu_h(-1) + epsilon_h;                      // Interest elasticity (markup) shocks
//mu_e=0;
//mu_h=0;

nu_h = rho_nuh*nu_h(-1) + epsilon_nu_h;                 // LTV shocks
nu_e = rho_nue*nu_e(-1) + epsilon_nu_e;                 // LTV shocks
//nu_h=0;
//nu_e=0;

xi_p = rho_p*xi_p(-1) + epsilon_p;                      // Cost-push shock
xi_t = rho_t*xi_t(-1) + epsilon_t;                      // capital-asset shock
xi_psi = rho_psi*xi_psi(-1) + epsilon_psi;              // Equity price shock
//xi_t=0;
//xi_psi=0;

//Measurement equations
y_obs = y - y(-1) + y_ss; // + me_y;
ygap_obs = ygap;
c_obs = c - c(-1) + c_ss; // + me_c;
v_obs = v - v(-1) + v_ss; // + me_v;
w_obs = w - w(-1) + w_ss; // + me_w;
h_obs = h + h_ss;         // + me_h;
pi_obs = pi + pi_ss;      // + me_p;

q_psi_obs = q_psi - q_psi(-1) + q_psi_ss; // + me_qpsi;
b_obs = b - b(-1) + b_ss;                 // + me_b;
l_h_obs = l_h - l_h(-1) + l_h_ss;         // + me_lh;
l_e_obs = l_e - l_e(-1) + l_e_ss;         // + me_le;
//k_BL_obs = k_BL - k_BL(-1) + k_BL_ss;
k_BL_obs = xi_t + k_BL_ss;

i_obs = i + i_ss;
//i_com_obs = i_com + i_com_ss;
//i_h_obs = i_h + i_h_ss;
//i_e_obs = i_e + i_e_ss;

//wage inflation
pi_wage = w - w(-1) + pi;                      // real wage definition // nominal wage inflation definition
winflationa = pi_wage*4;

//annualized inflation and interest rate
inflationa = pi*4;
interesta = i*4;

//output gap 
ygap = y - yn;

end;

// steady;
// check;


shocks;
var epsilon_z; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_g; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_v; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_i; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_p; stderr 0.005; // var(epsil) = 0.01^2
var epsilon_w; stderr 0.005; // var(epsil) = 0.01^2
var epsilon_b; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_psi; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_e; stderr 0.1; // var(epsil) = 0.01^2
var epsilon_h; stderr 0.1; // var(epsil) = 0.01^2
var epsilon_nu_h; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_nu_e; stderr 0.01; // var(epsil) = 0.01^2
var epsilon_t; stderr 0.01; // var(epsil) = 0.01^2
// var epsilon_k_BL; stderr 0.01; // var(epsil) = 0.01^2

%%% measurement errors % not used
%var me_y; stderr (0.0087/(100/10));        // (0.026752018/(100/15)) sqrt
%var me_c; stderr (0.0068/(100/10));        // (0.026752018/(100/15)) sqrt
%var me_v; stderr (0.006/(100/10));        // (0.026752018/(100/15)) sqrt
%var me_w; stderr (0.006/(100/10));        // (0.026752018/(100/15)) sqrt
%var me_h; stderr (0.006/(100/10));        // (0.026752018/(100/15)) sqrt
%var me_p; stderr (0.0025/(100/10));        // (0.026752018/(100/15)) sqrt
%var me_qpsi; stderr	(0.06/(100/15));     // (0.079364935/(100/25)) sqrt
%var me_b; stderr (0.016/(100/15));        // 0.007 (0.026752018/(100/15)) sqrt
%var me_lh; stderr (0.012/(100/10));        // 0.023613633 sqrt
%var me_le; stderr (0.012/(100/10));        // 0.035875444 sqrt

end;


estimated_params;

  //Households
gam, inv_gamma_pdf, 1.5, 0.5;           // coefficient of relative risk aversion
eta_h, gamma_pdf, 1.5, 0.5;             //labour disutility
phi_h, beta_pdf, 0.7, 0.1;              // habit formation
//phi_w, beta_pdf, 0.5, 0.1;            // 0.45 weight in borrowing constraint
Nu_h, beta_pdf, 0.5, 0.1;               // HH LTV ratio

  // Retailers and Unions
theta_R, beta_pdf, 0.7, 0.1;           // price stickiness
gamma_p, beta_pdf, 0.3, 0.1;            // degree of wage indexation
////varepsilon_p, inv_gamma_pdf, 7.667, inf;
theta_w, beta_pdf, 0.7, 0.1;           // wage stickiness
gamma_w, beta_pdf, 0.3, 0.1;           // degree of wage indexation

   //Entrepreneurs
gamma_e, inv_gamma_pdf, 1, 0.25;        // coefficient of relative risk aversion
//phi_k, beta_pdf, 0.5, 0.1;              // weight in borrowing constraint
Nu_e, beta_pdf, 0.5, 0.1;               // LTV ratio

   //Banks
////phi_psi, gamma_pdf, 0.3, 0.1;         // equity-to-capital ratio
Nu_B, beta_pdf, 0.2, 0.1;              // value-at-risk LTV ratio

   //Adjustment costs
kappa_k, gamma_pdf, 2, 1;              // leverage deviation cost
kappa_e, gamma_pdf, 5, 2;              // Entrepreneur loan rate adjustment costs
kappa_h, gamma_pdf, 5, 2;              // HH loan rate adjustment costs
kappa_v, gamma_pdf, 5, 0.5;              // physical capital adjustment costs

  //Monetary Policy function
////kappa_i, beta_pdf, 0.5, 0.1;         // Taylor rule coefficient on i
kappa_pi, gamma_pdf, 1.5, 0.1;        // Taylor rule coefficient on pi
kappa_y, beta_pdf, 0.5, 0.1;          // Taylor rule coefficient on y
////kappa_pich, beta_pdf, 0.35, 0.10;

   //Aggregate resource/flow of funds constraint
   
//K_eY, gamma_pdf, 10.7, 0.2;              //
//C_eY, inv_gamma_pdf, 0.04, 0.1;  
//CY, inv_gamma_pdf, 0.6734, 0.1;              
//K_BY, gamma_pdf, 0.1804, 0.02;              
//LY, gamma_pdf, 0.9, 0.3;
//PsiY, gamma_pdf, 0.7672, 0.2; 
//BY, gamma_pdf, 0.8, 0.05;
//HY, gamma_pdf, 0.61, 0.05;           // i.e. (1-alph)/X

//coefficients of AR(1) shocks

rho_z, beta_pdf, 0.5, 0.1;           // AR(1) parameter for productivity shock
rho_p, beta_pdf, 0.5, 0.1;           // AR(1) parameter for Cost-push shock
rho_w, beta_pdf, 0.5, 0.1;           // wage markup shock
rho_v, beta_pdf, 0.5, 0.1;           // investment shock
//rho_g, beta_pdf, 0.5, 0.1;         // exogenous spending shock
rho_i, beta_pdf, 0.5, 0.1;           // AR(1) parameter for MP shock

rho_b, beta_pdf, 0.5, 0.1;           // AR(1) parameter for deposit shock
rho_e, beta_pdf, 0.5, 0.1;           // AR(1) parameter for entrep shock
rho_h, beta_pdf, 0.5, 0.1;           // AR(1) parameter for HH interest elasticity  shock
rho_nuh, beta_pdf, 0.5, 0.1;
rho_nue, beta_pdf, 0.5, 0.1;         // AR(1) parameter for LTV shock, added ad hoc
rho_psi, beta_pdf, 0.75, 0.1;         // AR(1) parameter for Equity shock
rho_t, beta_pdf, 0.5, 0.1;           // AR(1) parameter for capital-asset shock shock

cgy, beta_pdf, 0.5, 0.1;             //weight on technology shock in exo spending process

//innovations of the 10 AR(1) shocks

stderr epsilon_z, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_g, inv_gamma_pdf, 0.01, inf;           // exogenous spending
stderr epsilon_v, inv_gamma_pdf, 0.01, inf;           // investment specific
stderr epsilon_i, inv_gamma_pdf, 0.01, inf;
stderr epsilon_p, inv_gamma_pdf, 0.01, inf;
stderr epsilon_w, inv_gamma_pdf, 0.01, inf;

stderr epsilon_b, inv_gamma_pdf, 0.01, inf;
stderr epsilon_psi, inv_gamma_pdf, 0.01, inf;
stderr epsilon_e, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_h, inv_gamma_pdf, 0.01, inf; 
stderr epsilon_nu_h, inv_gamma_pdf, 0.01, inf;
stderr epsilon_nu_e, inv_gamma_pdf, 0.01, inf;
stderr epsilon_t, inv_gamma_pdf, 0.01, inf;
// stderr epsilon_k_BL, inv_gamma_pdf, 0.01, inf;     // bank capital accumulation shock

end;

%stoch_simul(order=1, irf=20, graph_format=fig, nodisplay) y_obs inflationa interesta y c v h winflationa pi_wage;

// estimated_params_init(use_calibration);
// end;

 varobs y_obs c_obs v_obs w_obs h_obs pi_obs q_psi_obs b_obs l_h_obs l_e_obs i_obs k_BL_obs; //        i_com_obs   i_h_obs i_e_obs  
// varobs ygap_obs c_obs v_obs w_obs h_obs pi_obs q_psi_obs b_obs l_h_obs l_e_obs i_obs k_BL_obs; // i_com_obs   i_h_obs i_e_obs  


%% to run from scratch remove "mode_file= " and set "mode_compute=4"
estimation(datafile=HLEstimation12kbl, mh_nblocks=3, mh_replic=100000, mh_drop=0.5, mh_jscale=0.34,   // load_mh_file, 
           mode_compute=0, mode_file=HvL_HL16_rep13c_esttest1_mode, mode_check, bayesian_irf, irf=20, moments_varendo, conditional_variance_decomposition=[1:20],  // mode_file=HvL_HL16_rep_mode, 
           graph_format=fig) y_obs inflationa interesta ygap y c v h winflationa pi_wage i_com i_h i_e l_h l_e k_B k_BL b c_h c_e k_e q_k q_psi S_com S_h S_e l;

// Shock Decomposition option

shock_groups(name=groups);
'Demand' = epsilon_i, epsilon_v, epsilon_g;
'Supply' = epsilon_z, epsilon_p, epsilon_w;
'Credit demand' = epsilon_nu_h, epsilon_nu_e;
'Credit supply' = epsilon_h, epsilon_e, epsilon_t;  
'Savings' = epsilon_psi, epsilon_b;
//'Measurement' = epsilon_k_BL;         
end;

shock_groups(name=SWgroups);
'Technology' = epsilon_z;
'Markup' = epsilon_p, epsilon_w;
'Demand' = epsilon_v, epsilon_g;
'Monetary policy' = epsilon_i;
'Credit demand' = epsilon_nu_h, epsilon_nu_e;
'Credit supply' = epsilon_h, epsilon_e, epsilon_t;  
'Savings' = epsilon_psi, epsilon_b;
//'Measurement' = epsilon_k_BL;         
end;

%% shock_decomposition y_obs pi_obs S_com S_h S_e;
 shock_decomposition(colormap=jet) q_psi_obs l_h_obs l_e_obs ;  // b_obs
 plot_shock_decomposition(graph_format=(fig), use_shock_groups=groups, colormap=jet, steadystate) y_obs inflationa interesta ygap y yn S_com S_h S_e;

%% 

// //identification(parameter_set = prior_mode);
// identification(advanced=1, max_dim_cova_group=2, parameter_set = posterior_mean);