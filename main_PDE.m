% Accompanying code for the Paper:
% "Space--time discretization of the wave equation in a second-order-in-time
% formulation: a conforming, unconditionally stable method"
% Authors: M. Ferrari and I. Perugia
%
% This code computes an approximate solution of the wave problem:
% 
%     d^2_t U - Delta U = F,    x,t in Q_T = (0,L) x (0, T)
%     U(0,x) = dt U(0,x) = 0,   x   in (0,L)
%     U(t,0) = U(t,L) = 0,      t   in (0,T) 
%
% using the variational formulation:
% find U in L^2(0,T;H_0^1(0,L)) \cap H^1(0,T;L^2(0,L)) \cap H_{0,.}^2(0,T;H^{-1}(0,L))
%
%     (d_tt U, d_t W e^{-./T})_{L^2(Q_T)} 
%     + (d_t U(0,.), d_t W(0,.))_{L^2(Q_T)}
%     + (U, d_t W e^{-./T})_{L^2(Q_t)} 
%     = (F, d_t W e^{-./T})_{L^2(Q_T)}
%
% for all W in H^1(0,T;H_0^1(0,L)) \int H_{0,.}^2(0,T;L^2(0,L)).
%
% The code employs maximal regularity splines in space and arbitrary
% regularity splines in time for the numerical solution 

% and calculates errors in the norm L^2 and seminorms H^1, H^2 defined as
% || U ||_{L^2(Q_T)} := \int_0^T \int_0^L U^2(x,t) dx dt)^{1/2}, 
%  | U |_{H^1(Q_T)} := || dt U \|_{L^2(Q_T)}   + || dx U ||_{L^2(Q_T)} 
%  | U |_{H^2(Q_T)} := || d_tt U ||_{L^2(Q_T)} + || dx U ||_{L^2(Q_T)} 
%                                              + || dx dt U ||_{L^2(Q_T)}.

% Note: the code is far from being optimized

clc
clear
close all
format long

p = 4; % Degree of spline functions
T = 1; % Final time
L = 1; % Space interval [0,L]
r = 1; % Regularity splines

N_plot = 100; % Number of points for reconstructing the solution
Nq = 16; % Number of quadrature points
N1 = 2; % Starting resolution level (time and space intervals 2^(N1+1))
N2 = 3; % Ending resolution level (time and space intervals 2^(N2+1))
%also different mesh sizes for space and time can be employed

err_L2_app = zeros(1,N2-N1+1);
err_H1_app = zeros(1,N2-N1+1);
err_H2_app = zeros(1,N2-N1+1);

syms t x
%Analytical solution
%U = sin(pi*x)*sin(5/4*pi*t)^2; % U(0,t)=U(1,t)=0, U(x,0) = 0, dU(x,0)=0
%F = matlabFunction(diff(U,2,t) - diff(U,2,x));
%dtU = matlabFunction(diff(U,t));
%ddtU = matlabFunction(diff(U,2,t));
%dxU = matlabFunction(diff(U,x));
%ddxU = matlabFunction(diff(U,2,x));
%dxdtU = matlabFunction(diff(diff(U,1,x),1,t));
%U = matlabFunction(U);

%in good vectorial notation
U = @(t,x) sin(t*pi*5/4).^2*sin(x*pi);
F = @(t,x) pi^2*cos(t*pi*5/4).^2*sin(x.*pi)*(25/8)-pi^2*sin(t.*pi*5/4).^2*sin(x*pi)*17/8;
dtU =  @(t,x) pi*cos(t*pi*5/4).*sin(t*pi*5/4)*sin(x*pi)*5/2;
dxU = @(t,x) sin(t*pi*5/4).^2*pi*cos(x*pi);
ddtU = @(t,x) pi^2*cos(t*pi*5/4).^2*sin(x*pi)*25/8-pi^2*sin(t*pi*5/4).^2*sin(x*pi)*25/8;
ddxU = @(t,x) -pi^2*sin(t*pi*5/4).^2*sin(x*pi);
dtdxU = @(t,x) pi^2*cos(t*pi*5/4).*sin(t*pi*5/4)*5/2*cos(x*pi);

cont = 0;

for ii = N1:N2
    cont = cont + 1;
    N_x = 2^(ii+1);
    N_t = 2^(ii+1);
    
    knots_t = [zeros(r+1,1); repelem(linspace(0,T,N_t+1),p-r)'; ones(r+1,1)*T]';
    knots_x = [zeros(p,1); linspace(0,L,N_x+1)'; ones(p,1)*L]';

    t = ((0:N_t)/N_t)*T;
    x = ((0:N_x)/N_x)*L;
    h_t = t(2);
    h_x = x(2);

    B_t = mat_splines_exp(N_t,p,r,2,1,knots_t,T); % (dt^2 phi, dt phi e^{-./T})
    M_t = mat_splines_exp(N_t,p,r,0,1,knots_t,T); % (phi, dt phi e^{-./T})
    B_x = mat_splines(N_x,p,1,1,knots_x,L); % (dt phi, dt phi)
    M_x = mat_splines(N_x,p,0,0,knots_x,L); % (phi, phi)

    siz_t = N_t*(p-r)+r+1;
        F_vec = zeros((N_x+p)*siz_t,1);
        for j_t = 1 : siz_t
            for j_x = 1 : N_x+p
                for k_t = max(1,floor(j_t/(p-r))-(p+r)):min(j_t,N_t)
                    phi_t = zeros(Nq,1);
                    for k_x = max(1,j_x-p):min(j_x,N_x)
                        phi_x = zeros(Nq,1);
                        [xs_t,ps_t] = lgwt(Nq,t(k_t),t(k_t+1));
                        [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
                        for iii = 1:Nq
                            phi_t(iii) = sp_and_der(p,knots_t,j_t-1,xs_t(iii),1)*exp(-xs_t(iii)/T);
                            phi_x(iii) = sp_and_der(p,knots_x,j_x-1,xs_x(iii),0);
                        end
                        F_vec((j_t-1)*(N_x+p)+j_x) = F_vec((j_t-1)*(N_x+p)+j_x) + (phi_t.*ps_t)'*F(xs_t,xs_x')*(phi_x.*ps_x);
                    end
                end
            end
        end

    %removing initial condition in time and boundary in space
    cont_v = [];
    for j_t = 2 : siz_t
        for j_x = 2 : N_x+p-1
            cont_v = [cont_v (j_t-1)*(N_x+p)+j_x];
        end
    end
    F_vec = F_vec(cont_v);

    B_t = B_t(2:end,2:end); %zero initial condition
    M_t = M_t(2:end,2:end); %zero initial condition
    B_x = B_x(2:end-1,2:end-1); %zero boundary conditions
    M_x = M_x(2:end-1,2:end-1); %zero boundary conditions


    P_t = zeros(siz_t,siz_t);
    for jj_1 = 1:siz_t
        for jj_2 = 1:siz_t
            P_t(jj_1,jj_2) = sp_and_der(p,knots_t,jj_1-1,eps,1)*sp_and_der(p,knots_t,jj_2-1,eps,1);
        end
    end
    P_t = P_t(2:end,2:end);

    S = kron(sparse(B_t+P_t),sparse(M_x)) + kron(sparse(M_t),sparse(B_x));
    U_app_coeff = S\F_vec;
    
    t_plot = linspace(0.0001,T-0.0001,N_plot);
    x_plot = linspace(0.0001,L-0.0001,N_plot);

    U_ex = U(t_plot',x_plot);
    dtU_ex = dtU(t_plot',x_plot);
    dxU_ex = dxU(t_plot',x_plot);
    ddxU_ex = ddxU(t_plot',x_plot);
    ddtU_ex = ddtU(t_plot',x_plot);
    dtdxU_ex = dtdxU(t_plot',x_plot);
    U_app = zeros(size(U_ex));
    dtU_app = zeros(size(U_ex));
    dxU_app = zeros(size(U_ex));
    ddxU_app = zeros(size(U_ex));
    ddtU_app = zeros(size(U_ex));
    dtdxU_app = zeros(size(U_ex));

    %we add zero values at zero in time and at the boundary in space
    U_app_coeff_w_b = zeros((N_x+p)*siz_t,1);
    cont_v = 0;
    for j_t = 1 : siz_t
        for j_x = 1 : N_x+p
            if j_t == 1
                U_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = 0;
            elseif j_x == 1
                U_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = 0;
            elseif j_x == N_x+p
                U_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = 0;
            else
                cont_v = cont_v+1;
                U_app_coeff_w_b((j_t-1)*(N_x+p)+j_x) = U_app_coeff(cont_v);
            end
        end
    end
    U_app_coeff = U_app_coeff_w_b;

    for i_plot_t = 1 : N_plot
        for i_plot_x = 1 : N_plot
            for ind_t = 1 : siz_t
                if t_plot(i_plot_t) >= knots_t(ind_t) && t_plot(i_plot_t) < knots_t(ind_t+p+1)
                    for ind_x = 1 : N_x+p
                        if x_plot(i_plot_x) >= knots_x(ind_x) && x_plot(i_plot_x) < knots_x(ind_x+p+1)
                            U_app(i_plot_t,i_plot_x) = U_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                            dxU_app(i_plot_t,i_plot_x) = dxU_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),1);

                            dtU_app(i_plot_t,i_plot_x) = dtU_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                            ddtU_app(i_plot_t,i_plot_x) = ddtU_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),2)* ...
                                sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                            ddxU_app(i_plot_t,i_plot_x) = ddxU_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),2);

                            dtdxU_app(i_plot_t,i_plot_x) = dtdxU_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*(N_x+p)+ind_x)* ...
                                sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),1);
                        end
                    end
                end
            end
        end
    end

    % Errori in time
    err_L2_app(cont) = sqrt(sum(sum(abs(U_app-U_ex).^2))/(N_plot^2))/sqrt(sum(sum(abs(U_ex).^2))/(N_plot^2))
    err_H1_app(cont) = (sqrt(sum(sum(abs(dtU_app-dtU_ex).^2))/(N_plot^2))...
        + sqrt(sum(sum(abs(dxU_app-dxU_ex).^2))/(N_plot^2)))./...
        (sqrt(sum(sum(abs(dtU_ex).^2))/(N_plot^2)) + sqrt(sum(sum(abs(dxU_ex).^2))/(N_plot^2)));
    err_H2_app(cont) = (sqrt(sum(sum(abs(ddtU_app-ddtU_ex).^2))/(N_plot^2))...
        + sqrt(sum(sum(abs(ddxU_app-ddxU_ex).^2))/(N_plot^2))...
        + sqrt(sum(sum(abs(dtdxU_app-dtdxU_ex).^2))/(N_plot^2)))/...
        (sqrt(sum(sum(abs(ddtU_ex).^2))/(N_plot^2))+ sqrt(sum(sum(abs(ddxU_ex).^2))/(N_plot^2))...
        + sqrt(sum(sum(abs(dtdxU_ex).^2))/(N_plot^2)))

end

%error estimates
if N1 < N2 && N_x == N_t
    plot_error({}, err_L2_app,T./2.^(N1+1:N2+1),1);
    plot_error({}, err_H1_app,T./2.^(N1+1:N2+1),1);
    plot_error({'L2 error', 'H1 error','H2 error'}, err_H2_app, T./2.^(N1+1:N2+1),1);
    title(['Various errors for maximal regularity splines with p=' num2str(p)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
end
