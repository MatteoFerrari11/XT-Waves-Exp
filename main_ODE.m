% Accompanying code for the Paper:
% "Space--time discretization of the wave equation in a second-order-in-time
% formulation: a conforming, unconditionally stable method"
% Authors: M. Ferrari and I. Perugia
%
% This code computes an approximate solution of the ODE problem:
%
%     d^2_t u + mu * u = f,  t in (0, T)
%     u(0) = 0,  dt u(0) = 0
%
% using the variational formulation:
% find u in H^2(0,T) with u(0) = 0 such that
%
%     (d_tt u, d_t w e^{-./T})_{L^2(0,T)}
%     + d_t u(0) * d_t w(0)
%     + mu * (u, d_t w e^{-./T})_{L^2(0,T)}
%     = (f, d_t w e^{-./T})_{L^2(0,T)}
%
% for all w in H^2(0,T) with w(0) = 0.
%
% The code employs splines for the numerical solution
% and calculates relative errors in the following norms:
% L^2, H^1, H^2, L^inf, W^1_inf, W^2_inf

clc
clear
close all
format long

p = 4; % Degree of spline functions
r = 1; % Regularity of spline functions
T = 3; % Final time
mu = 1; % Parameter in the differential equation

N_plot = 10000; % Number of points for reconstructing the solution
Nq = 64; % Number of quadrature points
N1 = 2; % Starting resolution level (time intervals 2^(N1+1))
N2 = 6; % Ending resolution level (time intervals 2^(N2+1))

err_L2_app = zeros(1,N2-N1+1);
err_Linf_app = zeros(1,N2-N1+1);
err_W1inf_app = zeros(1,N2-N1+1);
err_W2ind_app = zeros(1,N2-N1+1);
err_H1_app = zeros(1,N2-N1+1);
err_H2_app = zeros(1,N2-N1+1);

syms t
u = exp(-t).*t.^2; % Analytical solution
du = matlabFunction(diff(u));
ddu = matlabFunction(diff(u,2));
f = matlabFunction(diff(u,2)+mu*u);
u = matlabFunction(u);

cont = 0;

for ii = N1:N2
    
    cont = cont+1;

    Nt = 2^(ii+1);
    knots = [zeros(r+1,1); repelem(linspace(0,T,Nt+1),p-r)'; ones(r+1,1)*T]';

    t = ((0:Nt)/Nt)*T;
    h = t(2);

    W1 = mat_splines_exp(Nt,p,r,2,1,knots,T); % (dt^2 phi, dt phi e^{-./T})
    W2 = mat_splines_exp(Nt,p,r,0,1,knots,T); % (phi, dt phi e^{-./T})

    siz = Nt*(p-r)+r+1;
    F = zeros(siz,1);
    for jj = 1:siz
        phi = zeros(Nq,1);
        for k = max([1,floor(jj/(p-r))-(p+r)]):min([jj,Nt])
            [xs,ps] = lgwt(Nq,t(k),t(k+1));
            for iii = 1:Nq
                phi(iii) = sp_and_der(p,knots,jj-1,xs(iii),1);
            end
            F(jj) = F(jj) + sum(f(xs).*phi.*ps.*exp(-xs/T));
        end
    end

    S_stab = zeros(siz,siz);
    for jj_1 = 1:siz
        for jj_2 = 1:siz
            S_stab(jj_1,jj_2) = sp_and_der(p,knots,jj_1-1,eps,1)*sp_and_der(p,knots,jj_2-1,eps,1);
        end
    end

    S = W1 + S_stab + mu*W2;

    S = S(2:end,2:end);  %zero initial conditions
    F = F(2:end); %zero initial conditions

    u_app_coeff = [0; S\F];

    x_plot = linspace(0.0001, T-0.0001, N_plot);
    u_ex = u(x_plot);
    du_ex = du(x_plot);
    ddu_ex = ddu(x_plot);
    u_app = zeros(size(u_ex));
    du_app = zeros(size(u_ex));
    ddu_app = zeros(size(u_ex));

    for i_plot = 1:N_plot
        for ind = 1:siz
            u_app(i_plot) = u_app(i_plot) + u_app_coeff(ind) * sp_and_der(p, knots,ind-1,x_plot(i_plot),0);
            du_app(i_plot) = du_app(i_plot) + u_app_coeff(ind) * sp_and_der(p,knots,ind-1,x_plot(i_plot),1);
            ddu_app(i_plot) = ddu_app(i_plot) + u_app_coeff(ind) * sp_and_der(p,knots,ind-1,x_plot(i_plot),2);
        end
    end

    if N1 == N2
        % Plot u_app and u_ex (solution)
        figure(1)
        plot(x_plot,real(u_app),'r','LineWidth',1.5)
        hold on
        plot(x_plot,real(u_ex),'b--','LineWidth',1.5)
        hold off
        title(['$u_h(t)$ and $u(t)$ for maximal regularity splines with p=' num2str(p)],'FontSize',14,'Interpreter','latex')
        xlabel('$t$','FontSize',12,'Interpreter','latex')
        legend({'$u_h(t)$', '$u(t)$'},'FontSize',16,'Location','best',...
            'Interpreter','latex')
        grid on

        % Plot du_app and du_ex (first derivative)
        figure(2)
        plot(x_plot,real(du_app),'r','LineWidth',1.5)
        hold on
        plot(x_plot,real(du_ex),'b--','LineWidth',1.5)
        hold off
        title(['$\partial_t u_h(t)$ and $\partial_t u(t)$ for maximal regularity splines with p=' num2str(p)],'FontSize',14,'Interpreter','latex')
        xlabel('$t$','FontSize',12,'Interpreter','latex')
        legend({'$\partial_t u_h(t)$','$\partial_t u(t)$'},'FontSize', 16,...
            'Location','best','Interpreter','latex')
        grid on

        % Plot the absolute error in u
        figure(3)
        loglog(x_plot,abs(real(u_app)-real(u_ex)),'r','LineWidth',1.5)
        title(['Absolute error $|u_h(t)-u(t)|$ for maximal regularity splines with p=' num2str(p)],'FontSize',14,'Interpreter','latex')
        xlabel('$t$','FontSize',12,'Interpreter','latex')
        ylabel('$|u_h(t)-u(t)|$','FontSize',12,'Interpreter','latex')
        grid on
    end

    err_L2_app(cont) = sqrt(T/length(x_plot)*sum(abs(u_app-u_ex).^2))/sqrt(T/length(x_plot)*sum(abs(u_ex).^2));
    err_Linf_app(cont) = max(abs(u_app-u_ex))/max(abs(u_ex));
    err_W1inf_app(cont) = max(abs(du_app-du_ex))/max(abs(du_ex));
    err_W2ind_app(cont) = max(abs(ddu_app-ddu_ex))/max(abs(ddu_ex));
    err_H1_app(cont) = sqrt(T/length(x_plot)*sum(abs(du_app-du_ex).^2))/sqrt(T/length(x_plot)*sum(abs(du_ex).^2));
    err_H2_app(cont) = sqrt(T/length(x_plot)*sum(abs(ddu_app-ddu_ex).^2))/sqrt(T/length(x_plot)*sum(abs(ddu_ex).^2));
end


if N1 < N2
    plot_error({}, err_L2_app,T./2.^(N1:N2),p-1);
    plot_error({}, err_H1_app,T./2.^(N1:N2),p-1);
    plot_error({}, err_H2_app, T./2.^(N1:N2),p-1);
    plot_error({}, err_Linf_app, T ./ 2.^(N1:N2),p-1);
    plot_error({}, err_W1inf_app, T ./ 2.^(N1:N2),p-1);
    plot_error({'L2 error', 'H1 error','H2 error','Linf error', 'W1inf error', 'W2inf error'}, err_W2ind_app, T ./ 2.^(N1:N2),p-1);
    title(['Various errors for regularity =' num2str(r) ' and splines with p=' num2str(p)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
end
