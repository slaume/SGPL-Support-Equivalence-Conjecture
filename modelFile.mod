option solver cplex;
option show_stats 1;
option cplex_options 'baropt bardisplay=1 barstart=1 comptol=1e-12 crossover=0'; 

#model parameters
param T;
param I;
param J;
param L {i in 1..I};
param A {i in 1..I, j in 0..J};
param c {i in 1..I};
param p {t in 1..T+1, j in 0..J};
param d {j in 0..J};
param f {j in 0..J};
param load {i in 1..I} := sum{t in 1..T, j in 0..J} p[t,j]*A[i,j]/c[i];
set RiU {i in 1..I, U in 0..2^(J+1)-1} 
	:= max{j in 0..J} (A[i,j]*floor((U mod 2^(j+1))/2^j))..L[i] 
		union {c[i]};

set Jset {i in 1..I} := union {j in 0..J} (if A[i,j]=1 then {j} else {});

#parameters for constructing random network instances
param choiceOfOD {od in 1..8};
param AallOD {i in 1..I, od in 1..8};
param delOD1;
param delOD2;
param delOD3;
param count;
param auxx;
param ODpairs {od in 1..8} symbolic;

#variables
var sigmaG {t in 1..T, i in 1..I, r in 1..L[i]+1} >= 0;
var zetaG {t in 1..T, i in 1..I, j in Jset[i], r in 1..L[i]+2} >= 0;
var muG {t in 1..T, j in 0..J} >= 0;
var phiNR {t in 1..T+1};
var VNR {t in 1..T+1, i in 1..I, k in 1..L[i]};
var alphaNR {t in 1..T+1, i in 1..I, U in 0..2^(J+1)-1};

#compact linear program (overline{P}_G)
maximize boundG: sum{t in 1..T, j in 1..J} p[t,j]*f[j]*muG[t,j];
subject to conG1 {i in 1..I, r in 1..L[i]-1}: 
	sigmaG[1,i,r] = 1;
subject to conG2 {t in 2..T, i in 1..I, r in 1..L[i]-1}:
	sigmaG[t,i,r] = sigmaG[t-1,i,r] 
	- sum{j in Jset[i]} p[t-1,j]*(zetaG[t-1,i,j,r]-zetaG[t-1,i,j,r+1]);
subject to conG3 {i in 1..I}:
	sigmaG[1,i,L[i]] + sigmaG[1,i,L[i]+1]*(c[i]-L[i]) = c[i]-L[i]+1;
subject to conG4 {t in 2..T, i in 1..I}:
	sigmaG[t,i,L[i]] + sigmaG[t,i,L[i]+1]*(c[i]-L[i])
	= sigmaG[t-1,i,L[i]] + sigmaG[t-1,i,L[i]+1]*(c[i]-L[i])
	- sum{j in Jset[i]} p[t-1,j]*zetaG[t-1,i,j,L[i]];
subject to conG5 {t in 1..T, i in 1..I}:
	sigmaG[t,i,1] <= 1;
subject to conG6 {t in 1..T, i in 1..I, j in Jset[i]}:
	muG[t,j] = zetaG[t,i,j,1];
subject to conG7 {t in 1..T, i in 1..I, j in Jset[i], r in 1..L[i]+1}:
	zetaG[t,i,j,r+1] <= zetaG[t,i,j,r];
subject to conG8 {t in 1..T, i in 1..I, j in Jset[i], r in 1..L[i]+1}:
	zetaG[t,i,j,r] <= sigmaG[t,i,r];
subject to conG9 {t in 1..T, i in 1..I, r in 1..L[i]}:
	sigmaG[t,i,r+1] <= sigmaG[t,i,r];
	
#partially reduced (tilde{D}_G)
minimize boundGdualNR: phiNR[1] 
	+ sum{i in 1..I, k in 1..L[i]-1} VNR[1,i,k]
	+ sum{i in 1..I} VNR[1,i,L[i]]*(c[i]-L[i]+1);
subject to conGdualNR1 {t in 1..T, U in 0..2^(J+1)-1}:
	sum{i in 1..I} alphaNR[t,i,U] >= phiNR[t+1]-phiNR[t] 
	+ sum{j in 0..J} p[t,j]*f[j]*floor((U mod 2^(j+1))/2^j);
subject to conGdualNR2 {t in 1..T, i in 1..I, U in 0..2^(J+1)-1, r in RiU[i,U]}:
	sum{k in 1..L[i]-1} (VNR[t,i,k]-VNR[t+1,i,k])*(if r>=k then 1 else 0)
	+ (VNR[t,i,L[i]]-VNR[t+1,i,L[i]])*max(0,r-L[i]+1)
	+ sum{j in 0..J} p[t,j]*floor((U mod 2^(j+1))/2^j)*
		(sum{k in 1..L[i]-1} VNR[t+1,i,k]*A[i,j]*(if r = k then 1 else 0)
			+ VNR[t+1,i,L[i]]*A[i,j]*(if r >= L[i] then 1 else 0))
	>= alphaNR[t,i,U];
subject to conGdualNR3 {i in 1..I, k in 1..L[i]}: 
	VNR[T+1,i,k] = 0;
subject to conGdualNR4:
	phiNR[T+1] = 0;
	
#problem definitions
problem G: boundG, sigmaG, zetaG, muG, 
	conG9, conG5, conG1, conG2, conG3, conG4, conG6, conG7, conG8;
problem GdualNR: boundGdualNR, VNR, alphaNR, phiNR, 
	conGdualNR1, conGdualNR2, conGdualNR3, conGdualNR4;
	



	
	
	
	

