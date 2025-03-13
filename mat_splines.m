function K = mat_splines(Nt,p,k1,k2,knots,T)

%input:
% Nt --> number of time intervals
%
% p --> degree of the maximal regularity splines
%
% (k1,k2) --> number of derivative in the matrix we want  to compute
% i.e. we have
% K[i,j] = (\partial_t^{k_1} \phi_j^p, \partial_t^{k_2} \phi_i^p)_{L^2(0,T)}
% with \phi_i^p and \phi_j^p basis function of the space of maximal
% regularity B-splines of degree p 
%
% knots --> set of knots for the definition of the basis functions
%
% T --> final time

%output: the matrix K \in R^{(Nt+p) \times (Nt+p)}  with entries
% K[i,j] = (dt^{k_1} \phi_j^p, dt^{k_2} \phi_i^p)_{L^2(0,T)}
% with \phi_i^p and \phi_j^p basis function of the space of maximal
% regularity B-splines of degree p 

t = ((0:Nt)/Nt)*T;
K = zeros(Nt+p,Nt+p);

Nq = 16; %for the quadrature rules

[x,w] = lgwt(Nq,-1,1);

for jj = 1 : Nt+p
    for ii = 1:Nt+p
        phi_1 = zeros(Nq,1);
        phi_2 = zeros(Nq,1);
        for k = max([1,jj-p,ii-p]):min([jj,ii,Nt])
            xs = (t(k+1)-t(k))/2*x + (t(k+1)+t(k))/2;
            ws = (t(k+1)-t(k))/2*w;
            for iii = 1 : Nq
                phi_1(iii) = sp_and_der(p,knots,ii-1,xs(iii),k1);
                phi_2(iii) = sp_and_der(p,knots,jj-1,xs(iii),k2);
            end
            K(jj,ii) = K(jj,ii) + sum(phi_2.*phi_1.*ws);
        end
    end
end
