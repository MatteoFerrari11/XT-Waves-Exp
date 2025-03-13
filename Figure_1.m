% Accompanying code for the Paper:
% "Space--time discretization of the wave equation in a second-order-in-time
% formulation: a conforming, unconditionally stable method"
% Authors: M. Ferrari and I. Perugia
%
% This code generates Figure 1

clc
clear
close all
format long

T = 5;
mu = 10^5;
N_plot = 10000;
Nq = 64;

err_H2_app = zeros(5,4);
err_H1_app = zeros(5,4);
err_L2_app = zeros(5,4);

syms t
u = exp(-t).*t.^2; % Analytical solution
du = matlabFunction(diff(u));
ddu = matlabFunction(diff(u,2));
f = matlabFunction(diff(u,2)+mu*u);
u = matlabFunction(u);

for p = 2:6
    r = 1; % Regularity of spline functions

    for ii = 2:5

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


        S = S(2:end,2:end);
        F = F(2:end);

        u_app_coeff = [0; S\F];

        x_plot = linspace(0.0001, T-0.0001, N_plot);
        u_ex = u(x_plot);
        du_ex = du(x_plot);
        ddu_ex = ddu(x_plot);
        u_app = zeros(size(du_ex));
        du_app = zeros(size(du_ex));
        ddu_app = zeros(size(du_ex));

        for i_plot = 1:N_plot
            for ind = 1:siz
                u_app(i_plot) = u_app(i_plot) + u_app_coeff(ind) * sp_and_der(p,knots,ind-1,x_plot(i_plot),0);
                du_app(i_plot) = du_app(i_plot) + u_app_coeff(ind) * sp_and_der(p,knots,ind-1,x_plot(i_plot),1);
                ddu_app(i_plot) = ddu_app(i_plot) + u_app_coeff(ind) * sp_and_der(p,knots,ind-1,x_plot(i_plot),2);
            end
        end

        err_H2_app(p-1,ii-1) = (inv(T)*sqrt(T/length(x_plot)*sum(abs(ddu_app-ddu_ex).^2))+sqrt(T/length(x_plot)*sum(abs(du_app-du_ex).^2)))/(inv(T)*sqrt(T/length(x_plot)*sum(abs(ddu_ex).^2))+sqrt(T/length(x_plot)*sum(abs(du_ex).^2)));
        err_L2_app(p-1,ii-1) = sqrt(T/length(x_plot)*sum(abs(u_app-u_ex).^2))/sqrt(T/length(x_plot)*sum(abs(u_ex).^2));
        err_H1_app(p-1,ii-1) = sqrt(T/length(x_plot)*sum(abs(du_app-du_ex).^2))/sqrt(T/length(x_plot)*sum(abs(du_ex).^2));

    end
    plot_error({},err_H2_app(p-1,:),T./2.^(3:6),1);
    plot_error({}, err_H1_app(p-1,:),T./2.^(3:6),2);
    plot_error({}, err_L2_app(p-1,:),T./2.^(3:6),3);
end


figure(1)
legend('p=2','p=3','p=4','p=5','p=6')
title('H2 error')
xlabel('h')

figure(2)
legend('p=2','p=3','p=4','p=5','p=6')
title('H1 error')
xlabel('h')

figure(3)
legend('p=2','p=3','p=4','p=5','p=6')
title('L2 error')
xlabel('h')