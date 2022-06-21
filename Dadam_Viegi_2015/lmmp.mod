% Monetary Policy In An Economy With High Structural Unemployment
% DSGE Model




var pi, ahat, uhat, chat, i, nhat, xhat, dm, m, l;
varexo ea, ed, em, el;


parameters alp, bet, lam, gam, e, B, M, x, g, u, delt, rhoa, rhod, rhom, rhol, psi, rho, net0, net1, h0, hL, hF, k0, kL, kF, phip, phic, phiu;


alp = 1; 						% a constant -- elasticity of hiring cost
bet = 0.99; 					% discount rate 
lam = 1/12; 					% price duration
gam = 0.5; 						% index of wage rigidity
e = 6; 							% elasticity of subsititution
B = 0.12; 						% level of hiring cost
M = e/(e-1); 						% optimal gross markup
x = 0.15; 						% labour market tightness index
g = B*x^alp; 					% hiring cost
u = 0.23; 						% unemployment
delt = u*x/((1-u)*(1-x)); 		% exogenous separation rate
psi = 1-(1-bet*(1-delt))*g*M; 	% Philips curve parameter
rhoa = 0.9; 					% persistence of the productivity process
rhom = 0.9;                     % persistence of the monetary shock process
rhod = 0.9;                     % persistence of the demand shock process
rhol = 0.9;                     % persistence of the employment shock process
rho = -log(bet);
net0 = (1-g*(1+alp))/(1-delt*g); 
net1 = g*(1-delt)*(1+alp*(1-x))/(1-delt*g);
h0 = (alp*g*M/delt)*(1+bet*(1-delt)*(1-delt)*(1-x))+bet*(1-delt)*g*M*(net1-net0);
hL = -(alp*g*M/delt)*(1-delt)*(1-x)-bet*(1-delt)*g*M*net1; 
hF = -bet*(1-delt)*g*M*((alp/delt)-net0);
k0 = lam*h0/(1-u); 
kL = -lam*hL/(1-u);
kF = lam*hF/(1-u); 
phip = 1.5; 
phic = 0.5; 
phiu = 0; 
 


model;


pi = bet*pi(+1)-k0*uhat+kL*uhat(-1)+kF*uhat(+1)-lam*psi*gam*ahat; 										
delt*xhat = nhat-(1-delt)*(1-x)*nhat(-1); 															
chat = ahat+((1-g)/(1-delt*g))*nhat+(g*(1-delt)/(1-delt*g))*nhat(-1)-(alp*g/(1-delt*g))*delt*xhat+l; 
chat = chat(+1)-(i-pi(+1)-rho)+dm; 																
uhat=-(1-u)*nhat; 																						
i = phip*pi+phic*chat+phiu*uhat+m;
ahat = rhoa*ahat(-1)+ea;
dm = rhod*dm(-1)+ed;
m = rhom*m(-1)+em;
l = rhol*l(-1)+el;


end;


varobs pi chat i nhat;

initval;
 
ahat=0; 
uhat=0; 
chat=0; 
i=0; 
nhat=0; 
xhat=0;
ea=0;
ed=0;
em=0;
el=0;
 
end;


estimated_params;

alp, beta_pdf, 0.9, 0.12;  					 
%lam, normal_pdf, 1/12, 0.5;
%e, normal_pdf, 6, 0.9; 							
B, beta_pdf, 0.2, 0.2; 	
x, beta_pdf, 0.7, 0.13;
gam, beta_pdf, 0.5, 0.25;

phip, normal_pdf, 1.5, 0.15; 
phic, normal_pdf, 0.125, 0.031; 
phiu, normal_pdf, 0, 0.021; 						

rhoa, beta_pdf, 0.8, 0.2;
rhod, beta_pdf, 0.8, 0.2;
rhom, beta_pdf, 0.8, 0.2;
rhol, beta_pdf, 0.8, 0.2;
 			 	 			
stderr ea, inv_gamma_pdf, 1, inf;
stderr ed, inv_gamma_pdf, 1, inf;
stderr em, inv_gamma_pdf, 1, inf;
stderr el, inv_gamma_pdf, 1, inf;

end;


estimation (datafile=he,nobs=79,first_obs=1,mh_replic=50000,mh_nblocks=5,
mh_drop=0.5,mh_jscale=3,mode_compute=6); 

shock_decomposition; 

