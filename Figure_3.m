% Accompanying code for the Paper:
% "Space--time discretization of the wave equation in a second-order-in-time
% formulation: a conforming, unconditionally stable method"
% Authors: M. Ferrari and I. Perugia
%
% This code generates Figure 3

clc
clear
close all
format long

T = 5;
N_plot = 10000;
Nq = 64;

err_H1_app = zeros(5,4,3);

syms t
u = sin(t).^2.*exp(-t).*t.^2; % Analytical solution
du = matlabFunction(diff(u));
u = matlabFunction(u);

for kk = 1:3
    for p = kk+1:6

        r = p-kk; % Regularity of spline functions

        for ii = 2:5

            Nt = 2^(ii+1);
            knots = [zeros(r+1,1); repelem(linspace(0,T,Nt+1),p-r)'; ones(r+1,1)*T]';

            t = ((0:Nt)/Nt)*T;
            h = t(2);

            Q = mat_splines_exp(Nt,p,r,1,2,knots,T); % (dt phi, dt phi e^{-./T})

            siz = Nt*(p-r)+r+1;
            F = zeros(siz,1);
            for jj = 1:siz
                phi = zeros(Nq,1);
                for k = max([1,floor(jj/(p-r))-(p+r)]):min([jj,Nt])
                    [xs,ps] = lgwt(Nq,t(k),t(k+1));
                    for iii = 1:Nq
                        phi(iii) = sp_and_der(p,knots,jj-1,xs(iii),2);
                    end
                    F(jj) = F(jj) + sum(du(xs).*phi.*ps.*exp(-xs/T));
                end
            end

            F = F(2:end);
            Q = Q(2:end,2:end);

            Last_Q = zeros(size(Q(end,:)));
            for jj = 1:siz
                Last_Q(jj) = sp_and_der(p,knots,jj-1,T-10*eps,1);
            end
            Last_Q = Last_Q(2:end);
            Q(end,:) = Last_Q;
            F(end) = du(T);

            u_app_coeff = [0; Q\F];

            x_plot = linspace(0.0001, T-0.0001, N_plot);
            du_ex = du(x_plot);
            du_app = zeros(size(du_ex));

            for i_plot = 1:N_plot
                for ind = 1:siz
                    du_app(i_plot) = du_app(i_plot) + u_app_coeff(ind) * sp_and_der(p,knots,ind-1,x_plot(i_plot),1);
                end
            end

            err_H1_app(p-kk,ii-1,kk) = sqrt(T/length(x_plot)*sum(abs(du_app-du_ex).^2))/sqrt(T/length(x_plot)*sum(abs(du_ex).^2))

        end
        plot_error({}, err_H1_app(p-kk,:,kk),T./2.^(3:6),kk);
    end
end

figure(1)
legend('p=2','p=3','p=4','p=5','p=6')
title('H1 error for spline (p,p-1)')
xlabel('h')

figure(2)
legend('p=3','p=4','p=5','p=6')
title('H1 error for spline (p,p-2)')
xlabel('h')

figure(3)
legend('p=4','p=5','p=6')
title('H1 error for spline (p,p-3)')
xlabel('h')