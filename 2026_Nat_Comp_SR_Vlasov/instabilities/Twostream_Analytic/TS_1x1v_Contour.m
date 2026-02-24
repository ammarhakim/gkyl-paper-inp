%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grant Johnson
%Started: 7/12/2022 (Redone with Newton, 12/15/2025)
%Calculate the contour of the warm Twostream (1x1v) with newton
    % Assumes c = 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Use diagnostics:
diagnostics = false;

% Wave numbers to solve for
Nk = 200;
k_vec = linspace(0.0,0.16,Nk);
omega_Im = zeros(Nk,1);

% Newton solve steps:
N_step = 15;

% parameters
T = 0.04; % vth = 0.2
n0 = 0.5; %*(1/1.063779); %Per beam due to shift
vb = 0.99;
q = -1.0;
m = 1.0;

% Plot the MJ/ compute density moment:
if (diagnostics) 
    plot_maxwellian(T,vb,m,n0);
    plot_contour_of_DR_1D(T,vb,q,m,n0);
end


% Newton solve for each k
for i = 1:Nk
    k = k_vec(i);
    omega_Im(i) = newton_solve(N_step,k,T,n0,vb,q,m);
end


% Plot the results of the newton solve:
plot(k_vec,omega_Im)
xlabel("k")
ylabel("\omega")

% Data for figure:
k_vec_final = k_vec(4:end)';
omega_final = omega_Im(4:end);






% Newton solver ( Use newton method to solve for the roots )
function omega_Im = newton_solve(N_step,k,T,n0,vb,q,m)

iter = 1;

% Find an inital guess
omega_Im0 = inital_omega_guess(k,T,vb,q,m,n0);

% Converge the root further
if omega_Im0 ~= 0 
    % Intial guess
    err = 1.0;
    
    % Main iteration step
    while (err > 1e-15 && iter < N_step) % && abs(omega_Im0) <= omega_Im0_max
    
    
        % Compute the function D(w_i) and D'(w_i) for newton
        disp = D_warm_rel(k,omega_Im0,T,vb,q,m,n0);
        disp_prime = D_warm_rel_prime(k,omega_Im0,T,vb,q,m,n0);
    
        % Compute the error
        err = abs(disp/ disp_prime);
        omega_Im_1 = omega_Im0 - disp/ disp_prime;
    
        % Print the progress:
        fprintf("(%d) (k = %1.1f), w_0 = %1.16e + i %1.16e, err: %1.16e\n", iter - 1, k, real(omega_Im0), imag(omega_Im0),err)
    
        % Set the old value to the new value
        omega_Im0 = omega_Im_1;
    
        % Increment 
        iter = iter + 1;
    
    end
end

% Set the solution:
omega_Im = abs(omega_Im0);

end


% Compute an inital guess of the imaginary part of the growth rate
function [omega_Im0] = inital_omega_guess(k,T,vb,q,m,n0)

% Compute for a range of omega the growth rate
omega_Im0_max = 0.07;
NV = 50;
dispersion_vec = zeros(1,NV);
omega_imag = linspace(0,omega_Im0_max,NV);

for i = 1:NV
    omega_Im_val  = omega_imag(i);
    dispersion_vec(i) = D_warm_rel(k,omega_Im_val,T,vb,q,m,n0);
end

% Loop through and find where the sign changes, choose the higher value
iter = 1;
root_found = 0;
omega_Im0 = 0;
dispersion_vec = real(dispersion_vec);
dispersion_vec(isnan(dispersion_vec)) = 0;
while iter < NV-1 && root_found == 0
    if ( dispersion_vec(iter) * dispersion_vec(iter+1) < 0 )
        root_found = 1;
        omega_Im0 = omega_imag(iter + 1);
    end
    iter = iter + 1;
end

end

%Calculate the DR.
function [D] = D_warm_rel(k,omega_im,T,vb,q,m0,n0)
gamma_b = 1/sqrt(1-vb^2);
K1 = besselk(1,m0/T);
C = n0 * gamma_b * q^2 / ( 2 * omega_im^2 * T * K1) ;

% Compute the integral in D(w_i)
format long
fun = @(p) ( p .* (((p./sqrt(1+p.^2)) + vb) ./ ((sqrt(1+p.^2)) + sqrt(-1) * k * p / omega_im )) ) .* exp((-gamma_b/T)*(sqrt(1+p.^2) + vb*p))+...
           ( p .* (((p./sqrt(1+p.^2)) - vb) ./ ((sqrt(1+p.^2)) + sqrt(-1) * k * p / omega_im )) ) .* exp((-gamma_b/T)*(sqrt(1+p.^2) - vb*p));
int_val = integral(fun,-inf,inf,'RelTol',1e-15,'AbsTol',1e-15);

D = 1 + C*int_val;
end


%Calculate the DR partial derivative.
function [D] = D_warm_rel_prime(k,omega_im,T,vb,q,m0,n0)
gamma_b = 1/sqrt(1-vb^2);
K1 = besselk(1,m0/T);
C = n0 * gamma_b * q^2 / ( 2 * K1 * T ) ;

% Compute the integral in D'(w_i)
format long
fun = @(p) ( 2 * omega_im + sqrt(-1) * k * (p./sqrt(1+p.^2)) ) .* ( p .* ((p./sqrt(1+p.^2)) + vb) ./ ( sqrt(1+p.^2) .* ( omega_im^2 + sqrt(-1) * k * ( p ./ sqrt(1+p.^2) ) * omega_im ).^2 ) ) .* exp((-gamma_b/T)*(sqrt(1+p.^2) + vb*p))+...
           ( 2 * omega_im + sqrt(-1) * k * (p./sqrt(1+p.^2)) ) .* ( p .* ((p./sqrt(1+p.^2)) - vb) ./ ( sqrt(1+p.^2) .* ( omega_im^2 + sqrt(-1) * k * ( p ./ sqrt(1+p.^2) ) * omega_im ).^2 ) ) .* exp((-gamma_b/T)*(sqrt(1+p.^2) - vb*p));
int_val = integral(fun,-inf,inf,'RelTol',1e-15,'AbsTol',1e-15);

D = - C*int_val;
end


% DIAGNOTIC 1: Calculate the MJ.
function [] = plot_maxwellian(T,vb,m0,n0)
gamma_b = 1/sqrt(1-vb^2);
K1 = besselk(1,m0/T);

% Function (MJ)
fun = @(p) n0 / ( 2 * m0 * K1) * exp((-gamma_b/T)*(sqrt(1+p.^2) - vb*p) );
density = integral(fun,-inf,inf,'RelTol',1e-15,'AbsTol',1e-15);

% Plotting function
NV = 1000;
p_vec = linspace(-15,15, NV);
f = subs(fun,p_vec);

plot(p_vec,f)
xlabel("u")
ylabel("f(u)");
fprintf("MJ Density (n0): %1.16e\n",density/gamma_b)
end


% DIAGNOSTIC 2: Compute the contour of the DR
function [] = plot_contour_of_DR_1D(T,vb,q,m,n0)

%Make the contour
k = 0.13;
omega_Im0_max = 0.07;
NV = 50;
dispersion_vec = zeros(1,NV);
omega_imag = linspace(0,omega_Im0_max,NV);

for i = 1:NV
    omega_Im0  = omega_imag(i);
    dispersion_vec(i) = D_warm_rel(k,omega_Im0,T,vb,q,m,n0);
end

figure();
plot(omega_imag,real(dispersion_vec))
xlabel("\gamma")
ylabel("D(\omega_I)")

end