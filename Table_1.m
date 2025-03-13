% Accompanying code for the Paper:
% "Space--time discretization of the wave equation in a second-order-in-time
% formulation: a conforming, unconditionally stable method"
% Authors: M. Ferrari and I. Perugia
%
% This code computes generates Table 1

clc
clear
close all
format short

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2;
T = 1;
mu = 10;

h_set_1 = zeros(8,1);
T_1_p_2_mu_10 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h_set_1(cont) = t(2);
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_1_p_2_mu_10(cont) = round(min(eig(S,H1)),3);
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2;
T = 1;
mu = 10^5;

T_1_p_2_mu_10_5 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_1_p_2_mu_10_5(cont) = round(min(eig(S,H1)),3);
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 3;
T = 1;
mu = 10;

T_1_p_3_mu_10 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_1_p_3_mu_10(cont) = round(min(eig(S,H1)),3);
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 3;
T = 1;
mu = 10^5;

T_1_p_3_mu_10_5 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_1_p_3_mu_10_5(cont) = round(min(eig(S,H1)),3);
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = round(h_set_1,3);
Tab_T_1 = table(h,T_1_p_2_mu_10,T_1_p_2_mu_10_5,T_1_p_3_mu_10,T_1_p_3_mu_10_5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2;
T = 3;
mu = 10;

h_set_3 = zeros(8,1);
T_3_p_2_mu_10 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h_set_3(cont) = t(2);
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_3_p_2_mu_10(cont) = round(min(eig(S,H1)),3);
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2;
T = 3;
mu = 10^5;

T_3_p_2_mu_10_5 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_3_p_2_mu_10_5(cont) = round(min(eig(S,H1)),3);
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 3;
T = 3;
mu = 10;

T_3_p_3_mu_10 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_3_p_3_mu_10(cont) = round(min(eig(S,H1)),3);
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 3;
T = 3;
mu = 10^5;

T_3_p_3_mu_10_5 = zeros(8,1);

cont = 0;

for ii = 2:9
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(p,1); linspace(0,T,Nt+1)'; ones(p,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);
    
    W1 = mat_splines_exp(Nt,p,p-1,2,1,knots,T);
    W2 = mat_splines_exp(Nt,p,p-1,0,1,knots,T);
    H1 = mat_splines(Nt,p,1,1,knots,T); % (dt phi, dt phi)

    S = W1 + mu*W2;
    S = S(2:end,2:end); %zero initial conditions
    S(1,1) = S(1,1) + p^2/h^2;

    H1 = H1(2:end,2:end);
    S = (S+S')/2;
    T_3_p_3_mu_10_5(cont) = round(min(eig(S,H1)),3);
  
end

h = round(h_set_3,3);
Tab_T_3 = table(h,T_3_p_2_mu_10,T_3_p_2_mu_10_5,T_3_p_3_mu_10,T_3_p_3_mu_10_5)
