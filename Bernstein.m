clear all;
close all;

% number of modes
n=20;
% size of the region (i.e. scaling of a region [0,1]x[0,1])
alpha=3;
alpha1=alpha;
alpha2=alpha;
% position of the fixed point in the region
% 0: left border;0.5: middle; 1: right border 
gamma_ref1=0.5;
% 0: bottom border;0.5: middle; 1: top border 
gamma_ref2=0.5;


% rotating dynamics x1dot=-x1-x1^2*x2-x2^3; x2dot=-x2+x1*x2^2+x1^3 
% (equivalent to dotr=-r; dottheta=r^2)
pt_fix1=0;
pt_fix2=0;
J=[-1 0;0 -1];
% coefficients [B00 B01 B02 B03 B10 B11 B12 B13 B20 B21 B22 B23 B30 B31 B32 B33]
pol_F1=[0 0 0 -1 -1 0 0 0 0 -1 0 0 0 0 0 0]';
pol_F2=[0 -1 0 0 0 0 1 0 0 0 0 0 1 0 0 0]';


%% scaling

[u,v]=eig(J');
lambda1=v(2,2);
lambda2=v(1,1);
vec1=u(:,2);
vec2=u(:,1);

pt_ref1=pt_fix1;
pt_ref2=pt_fix2;

m=sqrt(length(pol_F1))-1; 

pol_F1b=zeros(length(pol_F1),1);
pol_F2b=zeros(length(pol_F1),1);

for k=1:length(pol_F1)
    k1=ceil(k/(m+1))-1;
    k2=mod(k-1,m+1);
    for i=0:k1
        for j=0:k2
            pol_F1b((k1-i)*(m+1)+(k2-j)+1)=pol_F1b((k1-i)*(m+1)+(k2-j)+1)+pol_F1(k)*alpha1^k1*alpha2^k2*nchoosek(k1,i)*nchoosek(k2,j)*(-(gamma_ref1-pt_ref1))^(i+j);
            pol_F2b((k1-i)*(m+1)+(k2-j)+1)=pol_F2b((k1-i)*(m+1)+(k2-j)+1)+pol_F2(k)*alpha1^k1*alpha2^k2*nchoosek(k1,i)*nchoosek(k2,j)*(-(gamma_ref2-pt_ref2))^(i+j);
        end
    end
end

pol_F1b=pol_F1b/alpha;
pol_F2b=pol_F2b/alpha;

%% dynamics in Bernstein polynomials

U1=zeros(m+1,m+1);
for i=0:m
    for j=0:i
        U1(i+1,j+1)=nchoosek(m,i)*nchoosek(i,j)*(-1)^(i-j);
    end
end
U=kron(U1,U1);
b_F1=inv(U)*pol_F1b;
b_F2=inv(U)*pol_F2b;

%% computation of the matrices

% matrices for derivatives
D=[-n*eye(n) zeros(n,1)]+[zeros(n,1) n*eye(n)];
D1=kron(D,eye(n+1));
D2=kron(eye(n+1),D);

% raising the polynomial degree for the right hand side (from n*n to (m+n)*(m+n))
T=zeros(m+n+1,n+1);
for j=0:n
    for i=j:j+m
    T(i+1,j+1)=nchoosek(n,j)*nchoosek(m,i-j)/nchoosek(m+n,i);
    end
end
T_tot=kron(T,T);

% raising the polynomial degree from (n-1)*n to n*n
T=zeros(n+1,n);
for j=0:n-1
    for i=j:j+1
    T(i+1,j+1)=nchoosek(n-1,j)*nchoosek(1,i-j)/nchoosek(n,i);
    end
end
T1=kron(T,eye(n+1));
T2=kron(eye(n+1),T);

% matrices for multplication
M1=zeros((m+n+1)*(m+n+1),(n+1)*(n+1));
M2=zeros((m+n+1)*(m+n+1),(n+1)*(n+1));
for k1=0:m
    for k2=0:m
        M1k=zeros(m+n+1,n+1);
        M2k=zeros(m+n+1,n+1);
        for i=1:m+n+1
            for j=1:n+1
                if j==i-k1
                    M1k(i,j)=nchoosek(n,j-1)*nchoosek(m,k1)/nchoosek(n+m,i-1);
                end
                if j==i-k2
                    M2k(i,j)=nchoosek(n,j-1)*nchoosek(m,k2)/nchoosek(n+m,i-1);
                end
            end
        end
        M1=M1+b_F1(k1*(m+1)+k2+1)*kron(M1k,M2k);
        M2=M2+b_F2(k1*(m+1)+k2+1)*kron(M1k,M2k);
    end
end

%%
% additional constraint: zero at the fixed point
C1=zeros(1,(n+1)^2);
for k1=0:n
    for k2=0:n
        C1(k1*(n+1)+k2+1)=nchoosek(n,k1)*gamma_ref1^k1*(1-gamma_ref1)^(n-k1)*nchoosek(n,k2)*gamma_ref2^k2*(1-gamma_ref2)^(n-k2);
    end
end

% additional constraints: gradient at the fixed point is equal to left eigenvector of Jacobian matrix
% derivative with respect to x1
C2=zeros(1,(n+1)^2);
for k1=1:n-1
    for k2=0:n
            C2(k1*(n+1)+k2+1)=(nchoosek(n,k1)*(k1*gamma_ref1^(k1-1)*(1-gamma_ref1)^(n-k1)-gamma_ref1^k1*(n-k1)*(1-gamma_ref1)^(n-k1-1)))*nchoosek(n,k2)*gamma_ref2^k2*(1-gamma_ref2)^(n-k2);
    end
end
k1=0;
for k2=0:n
    C2(k1*(n+1)+k2+1)=-n*(1-gamma_ref1)^(n-1)*nchoosek(n,k2)*gamma_ref2^k2*(1-gamma_ref2)^(n-k2);
end
k1=n;
for k2=0:n
    C2(k1*(n+1)+k2+1)=n*gamma_ref1^(n-1)*nchoosek(n,k2)*gamma_ref2^k2*(1-gamma_ref2)^(n-k2);
end

% derivative with respect to x2
C3=zeros(1,(n+1)^2);
for k1=0:n
    for k2=1:n-1
            C3(k1*(n+1)+k2+1)=nchoosek(n,k1)*gamma_ref1^k1*(1-gamma_ref1)^(n-k1)*(nchoosek(n,k2)*(k2*gamma_ref2^(k2-1)*(1-gamma_ref2)^(n-k2)-gamma_ref2^k2*(n-k2)*(1-gamma_ref2)^(n-k2-1)));
    end
end
k2=0;
for k1=0:n
    C3(k1*(n+1)+k2+1)=nchoosek(n,k1)*gamma_ref1^k1*(1-gamma_ref1)^(n-k1)*(-n*(1-gamma_ref2)^(n-1));
end
k2=n;
for k1=0:n
    C3(k1*(n+1)+k2+1)=nchoosek(n,k1)*gamma_ref1^k1*(1-gamma_ref1)^(n-k1)*(n*gamma_ref2^(n-1));
end

%% computation of the eigenfunctions

n_weight=1;

A1=[M1*T1*D1+M2*T2*D2-lambda1*T_tot;C1;n_weight*C2;n_weight*C3];
b1=[zeros((m+n+1)*(m+n+1),1);0;n_weight*vec1];
A2=[M1*T1*D1+M2*T2*D2-lambda2*T_tot;C1;n_weight*C2;n_weight*C3];
b2=[zeros((m+n+1)*(m+n+1),1);0;n_weight*vec2];

b_phi1=pinv(A1)*b1;
b_phi2=pinv(A2)*b2;

test1=norm(A1*b_phi1-b1)/norm(b1);
test2=max(abs(A1*b_phi1-b1))/max(abs(b1));

%% reconstruction

x_interv=linspace(0,1,100);
y_interv=linspace(0,1,100);

[x0,y0]=meshgrid(x_interv,y_interv);

phi1=zeros(size(x0));
phi2=zeros(size(x0));

for k1=0:n
    for k2=0:n
        phi1=phi1+b_phi1(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
        phi2=phi2+b_phi2(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
    end
end

% grad_phi1_x1=zeros(size(x0));
% grad_phi2_x1=zeros(size(x0));
% for k1=1:n-1
%     for k2=0:n
%             grad_phi1_x1=grad_phi1_x1+b_phi1(k1*(n+1)+k2+1).*(nchoosek(n,k1).*(k1.*x0.^(k1-1).*(1-x0).^(n-k1)-x0.^k1.*(n-k1).*(1-x0).^(n-k1-1))).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
%             grad_phi2_x1=grad_phi2_x1+b_phi2(k1*(n+1)+k2+1).*(nchoosek(n,k1).*(k1.*x0.^(k1-1).*(1-x0).^(n-k1)-x0.^k1.*(n-k1).*(1-x0).^(n-k1-1))).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
%     end
% end
% k1=0;
% for k2=0:n
%     grad_phi1_x1=grad_phi1_x1-b_phi1(k1*(n+1)+k2+1).*n.*(1-x0).^(n-1).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
%     grad_phi2_x1=grad_phi2_x1-b_phi2(k1*(n+1)+k2+1).*n.*(1-x0).^(n-1).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
% end
% k1=n;
% for k2=0:n
%     grad_phi1_x1=grad_phi1_x1+b_phi1(k1*(n+1)+k2+1).*n.*x0.^(n-1).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
%     grad_phi2_x1=grad_phi2_x1+b_phi2(k1*(n+1)+k2+1).*n.*x0.^(n-1).*nchoosek(n,k2).*y0.^k2.*(1-y0).^(n-k2);
% end
% 
% grad_phi1_x2=zeros(size(x0));
% grad_phi2_x2=zeros(size(x0));
% for k1=0:n
%     for k2=1:n-1
%             grad_phi1_x2=grad_phi1_x2+b_phi1(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*(nchoosek(n,k2).*(k2.*y0.^(k2-1).*(1-y0).^(n-k2)-y0.^k2.*(n-k2).*(1-y0).^(n-k2-1)));
%             grad_phi2_x2=grad_phi2_x2+b_phi2(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*(nchoosek(n,k2).*(k2.*y0.^(k2-1).*(1-y0).^(n-k2)-y0.^k2.*(n-k2).*(1-y0).^(n-k2-1)));
%     end
% end
% k2=0;
% for k1=0:n
%     grad_phi1_x2=grad_phi1_x2+b_phi1(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*(-n.*(1-y0).^(n-1));
%     grad_phi2_x2=grad_phi2_x2+b_phi2(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*(-n.*(1-y0).^(n-1));
% end
% k2=n;
% for k1=0:n
%     grad_phi1_x2=grad_phi1_x2+b_phi1(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*(n.*y0.^(n-1));
%     grad_phi2_x2=grad_phi2_x2+b_phi2(k1*(n+1)+k2+1).*nchoosek(n,k1).*x0.^k1.*(1-x0).^(n-k1).*(n.*y0.^(n-1));
% end
% 
% 
% F1_x=zeros(size(x0));
% F2_x=zeros(size(x0));
% for k1=0:m
%     for k2=0:m
%         F1_x=F1_x+pol_F1b(k1*(m+1)+k2+1).*x0.^k1.*y0.^k2;
%         F2_x=F2_x+pol_F2b(k1*(m+1)+k2+1).*x0.^k1.*y0.^k2;
%     end
% end
% 
% phi1_modif=phi1-C1*b_phi1;
% phi2_modif=phi2-C1*b_phi2;
% % grad_V_x1=real(phi1_modif).*real(grad_phi1_x1)+imag(phi1_modif).*imag(grad_phi1_x1)+real(phi2_modif).*real(grad_phi2_x1)+imag(phi2_modif).*imag(grad_phi2_x1);
% % grad_V_x2=real(phi1_modif).*real(grad_phi1_x2)+imag(phi1_modif).*imag(grad_phi1_x2)+real(phi2_modif).*real(grad_phi2_x2)+imag(phi2_modif).*imag(grad_phi2_x2);
% grad_V_x1=sign(phi1_modif).*grad_phi1_x1+sign(phi2_modif).*grad_phi2_x1;
% grad_V_x2=sign(phi1_modif).*grad_phi1_x2+sign(phi2_modif).*grad_phi2_x2;
% Vdot=grad_V_x1.*F1_x+grad_V_x2.*F2_x;
% test=Vdot<0;
% test=+test;

%% plot

figure()
hold on
contour((x0-(gamma_ref1-pt_ref1))*alpha1,(y0-(gamma_ref2-pt_ref2))*alpha2,abs(phi1),20)
box on
plot(0,0,'ok','Linewidth',5,'MarkerSize',5)
xlabel('$x_1$','Interpreter','latex','Fontsize',20,'Rotation',0)
ylabel('$x_2$','Interpreter','latex','Fontsize',20,'Rotation',0)
colorbar

figure()
hold on
contour((x0-(gamma_ref1-pt_ref1))*alpha1,(y0-(gamma_ref2-pt_ref2))*alpha2,abs(phi2),20)
box on
plot(0,0,'ok','Linewidth',5,'MarkerSize',5)
xlabel('$x_1$','Interpreter','latex','Fontsize',20,'Rotation',0)
ylabel('$x_2$','Interpreter','latex','Fontsize',20,'Rotation',0)
colorbar

% figure()
% pcolor((x0-(gamma_ref1-pt_ref1))*alpha1,(y0-(gamma_ref2-pt_ref2))*alpha2,test)
% hold on
% contour((x0-(gamma_ref1-pt_ref1))*alpha1,(y0-(gamma_ref2-pt_ref2))*alpha2,abs(phi1)+abs(phi2),[0:0.1:1.],'k')
