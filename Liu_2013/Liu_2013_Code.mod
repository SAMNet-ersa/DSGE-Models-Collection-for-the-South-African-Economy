//code for liu4_Price puzzle paper (estimation)
//code for the model (no cost channel) taken from liu1_US money (model 2: interest rate rule only)


//AR(1) monetary policy shock
var y c i pie r n w r_k k p_k mc xi_z xi_x x xi_i; 

//IID monetary policy shock
//var y c i pie r n w r_k k p_k mc xi_z xi_x x;

varexo epsilon_z epsilon_x epsilon_p epsilon_w epsilon_i; 


parameters delta ky h eta_c psi alpha beta nu gamma_p theta_p gamma_w theta_w lambda_w eta_n kappa_i kappa_pie kappa_y rho_i rho_z rho_x gamma;




//parameters fixed prior estimation
delta= 0.019; 
ky = 1.7; 
alpha = 0.26; 
beta = 0.99; 
lambda_w = 21; //varphi_w in the paper


model(linear);
k = (1-delta)*k(-1) + delta*x; 
y = (1-delta*ky)*c + delta*ky*x; 
c = (h/(1+h))*c(-1) + (1/(1+h))*c(+1) - ((1-h)/((1+h)*eta_c))*(r);
i = r + pie(+1);
//n = -w + (1+psi)*r_k + k(-1); 
n = -(w-gamma*i) + (1+psi)*r_k + k(-1); //with cost channel, then you need to add gamma into your parameter
y = xi_z + alpha*(psi*r_k + k(-1))+(1-alpha)*n; 
x = (1/(1+beta))*x(-1)+(beta/(1+beta))*x(+1)+(1/(1+beta)*nu)*p_k+xi_x; 
p_k = pie(+1)-r+(1-beta*(1-delta))*r_k(+1)+beta*(1-delta)*p_k(+1); 
pie = (1/(1+beta*gamma_p))*(gamma_p*pie(-1)+beta*pie(+1)+((1-beta*theta_p)*(1-theta_p))/(theta_p)*mc)+epsilon_p; 
//mc = alpha*r_k+(1-alpha)*w-xi_z; 
mc = alpha*r_k+(1-alpha)*(w+gamma*i)-xi_z; //with cost channel, need to add gamma into your parameter
w = (1/(1+beta))*( w(-1)+beta*w(+1)+gamma_w*(pie(-1))-(1+beta*gamma_w)*pie+beta*(pie(+1))-(((1-beta*theta_w)*(1-theta_w))/((1+lambda_w*eta_n)*theta_w))*(w-eta_n*n-(eta_c/(1-h))*(c-h*(c(-1))))) + epsilon_w; 


//AR(1) monetary policy shock
i = kappa_i*i(-1) + (1-kappa_i)*(kappa_pie*pie + kappa_y*y)+ xi_i; // consider pie(+1)

//IID monetary policy shock
//i = kappa_i*i(-1) + (1-kappa_i)*(kappa_pie*pie+kappa_y*y)+ epsilon_i; 


xi_i = rho_i*xi_i(-1) + epsilon_i; 
xi_x = rho_x*xi_x(-1) + epsilon_x; 
xi_z = rho_z*xi_z(-1) + epsilon_z; 

end;




estimated_params;

eta_c, gamma_pdf, 2.1, 0.1;
eta_n, normal_pdf, 1.5, 0.75; 
h, beta_pdf, 0.7, 0.05; 
psi, gamma_pdf, 20, 4.14; 
nu, gamma_pdf, 2, 0.4; 

//gamma = 1; // to conduct a comparison study set gamma= 1, 0.7, 0.2, 0. experience show pie increases as long as gamma is greater than 3.6.
gamma, beta_pdf, 0.7, 0.1;

//given values of gamma and psi are very critical guarantee the BK condition is satisfied.

theta_p, beta_pdf, 0.50, 0.10; 
theta_w, beta_pdf, 0.75, 0.05; 
gamma_p, beta_pdf, 0.55, 0.15; 
gamma_w, beta_pdf, 0.5, 0.20; 

// coefficients of AR(1) shocks; NB: price and wage mark-up shocks are assumed to be serially uncorrelated   
rho_x,  beta_pdf, 0.75, 0.10; 
rho_i,  beta_pdf, 0.55, 0.05; 
rho_z,  beta_pdf, 0.95, 0.02; 

//coefficients in monetary policy function 

kappa_i, normal_pdf, 0.60, 0.05;
kappa_pie, normal_pdf, 1.8, 0.05;
kappa_y, normal_pdf, 0.25, 0.02;


//innovations of the shocks [you need to play around with the mean (in particular) and standard error]
stderr epsilon_i,  inv_gamma_pdf, 0.10, inf;
stderr epsilon_z,  inv_gamma_pdf, 0.10, inf;
stderr epsilon_x,  inv_gamma_pdf, 0.10, inf;
stderr epsilon_p,  inv_gamma_pdf, 0.10, inf;
stderr epsilon_w,  inv_gamma_pdf, 0.10, inf;

end;


varobs y n w pie i;  
estimation(datafile=Data_liu4e, mh_nblocks=2, mh_replic=20000, mh_jscale=0.5, mh_drop=0.5, prior_trunc=1e-100, mode_check, bayesian_irf, irf=20) y n w pie i;
stoch_simul(order=1, ar=12);