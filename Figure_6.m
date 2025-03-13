% Accompanying code for the Paper:
% "Space--time discretization of the wave equation in a second-order-in-time
% formulation: a conforming, unconditionally stable method"
% Authors: M. Ferrari and I. Perugia
%
% This code generates Figure 6. It is slow

clc
clear
close all
format long

T = 1; % Final time
L = 1; % Space interval [0,L]
N_plot = 100; % Number of points for reconstructing the solution
Nq = 16; % Number of quadrature points
N1 = 2; % Starting resolution level (time and space intervals 2^(N1+1))
N2 = 5; % Ending resolution level (time and space intervals 2^(N2+1))
%also different mesh sizes for space and time can be employed

syms t x
%Analytical solution
U = @(t,x) sin(t*pi*5/4).^2*sin(x*pi);
F = @(t,x) pi^2*cos(t*pi*5/4).^2*sin(x.*pi)*(25/8)-pi^2*sin(t.*pi*5/4).^2*sin(x*pi)*17/8;
dtU =  @(t,x) pi*cos(t*pi*5/4).*sin(t*pi*5/4)*sin(x*pi)*5/2;
dxU = @(t,x) sin(t*pi*5/4).^2*pi*cos(x*pi);
ddtU = @(t,x) pi^2*cos(t*pi*5/4).^2*sin(x*pi)*25/8-pi^2*sin(t*pi*5/4).^2*sin(x*pi)*25/8;
ddxU = @(t,x) -pi^2*sin(t*pi*5/4).^2*sin(x*pi);
dtdxU = @(t,x) pi^2*cos(t*pi*5/4).*sin(t*pi*5/4)*5/2*cos(x*pi);

for p = 2:3 %Degree of spline functions
    
    err_L2_app = zeros(1,N2-N1+1);
    err_H1_app = zeros(1,N2-N1+1);
    err_H2_app = zeros(1,N2-N1+1);

    cont = 0;

    for ii = N1:N2

        cont = cont + 1;
        N_x = 2^(ii+1);
        N_t = 2^(ii+1);

        knots_t = [zeros(p,1); linspace(0,T,N_t+1)'; ones(p,1)*T]';
        knots_x = [zeros(p,1); linspace(0,L,N_x+1)'; ones(p,1)*L]';

        t = ((0:N_t)/N_t)*T;
        x = ((0:N_x)/N_x)*L;
        h_t = t(2)
        h_x = x(2);

        B_t = mat_splines_exp(N_t,p,p-1,2,1,knots_t,T); % (dt^2 phi, dt phi e^{-./T})
        M_t = mat_splines_exp(N_t,p,p-1,0,1,knots_t,T); % (phi, dt phi e^{-./T})
        B_x = mat_splines(N_x,p,1,1,knots_x,L); % (dt phi, dt phi)
        M_x = mat_splines(N_x,p,0,0,knots_x,L); % (phi, phi)

        F_vec = zeros((N_x+p)*(N_t+p),1);
        for j_t = 1 : N_t+p
            for j_x = 1 : N_x+p
                for k_t = max(1,j_t-p):min(j_t,N_t)
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
        for j_t = 2 : N_t+p
            for j_x = 2 : N_x+p-1
                cont_v = [cont_v (j_t-1)*(N_x+p)+j_x];
            end
        end
        F_vec = F_vec(cont_v);

        B_t = B_t(2:end,2:end); %zero initial condition
        M_t = M_t(2:end,2:end); %zero initial condition
        B_x = B_x(2:end-1,2:end-1); %zero boundary conditions
        M_x = M_x(2:end-1,2:end-1); %zero boundary conditions

        P_t = zeros(size(B_t));
        P_t(1,1) = p^2/h_t^2;

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
        U_app_coeff_w_b = zeros((N_x+p)*(N_t+p),1);
        cont_v = 0;
        for j_t = 1 : N_t+p
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
                for ind_t = 1 : N_t+p
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

        err_L2_app(cont) = sqrt(sum(sum(abs(U_app-U_ex).^2))/(N_plot^2))/sqrt(sum(sum(abs(U_ex).^2))/(N_plot^2))
        err_H1_app(cont) = (sqrt(sum(sum(abs(dtU_app-dtU_ex).^2))/(N_plot^2))...
            + sqrt(sum(sum(abs(dxU_app-dxU_ex).^2))/(N_plot^2)))./...
            (sqrt(sum(sum(abs(dtU_ex).^2))/(N_plot^2)) + sqrt(sum(sum(abs(dxU_ex).^2))/(N_plot^2)))
        err_H2_app(cont) = (sqrt(sum(sum(abs(ddtU_app-ddtU_ex).^2))/(N_plot^2))...
            + sqrt(sum(sum(abs(ddxU_app-ddxU_ex).^2))/(N_plot^2))...
            + sqrt(sum(sum(abs(dtdxU_app-dtdxU_ex).^2))/(N_plot^2)))/...
            (sqrt(sum(sum(abs(ddtU_ex).^2))/(N_plot^2))+ sqrt(sum(sum(abs(ddxU_ex).^2))/(N_plot^2))...
            + sqrt(sum(sum(abs(dtdxU_ex).^2))/(N_plot^2)))

    end

    plot_error({}, err_L2_app,T./2.^(N1+1:N2+1),p-1);
    plot_error({}, err_H1_app,T./2.^(N1+1:N2+1),p-1);
    plot_error({'L2 error', 'H1 error','H2 error'}, err_H2_app, T./2.^(N1+1:N2+1),p-1);
    title(['Various errors for maximal regularity splines with p=' num2str(p)],'FontSize',16,'Interpreter','latex')
    xlabel('h')

end