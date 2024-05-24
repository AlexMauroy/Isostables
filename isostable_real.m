% This code computes the isostables of the FitzHugh-Nagumo model (case of real eigenvalues), using forward integration. It runs in about 2 minutes 
% on a standard computer.
% The code can be easily adapted to other models with a stable fixed point. The method also works with models of higher dimensions.

% The method is based on the results presented in "A. Mauroy, I. Mezic, and J. Moehlis, Isostables, isochrons, and Koopman spectrum 
% for the action-angle reduction of stable fixed point dynamics, Physica D, vol. 261, pp. 19-30, 2013"

% For more information, please email me at "alexandre.mauroy 'at' unamur.be"

% Written by Alexandre Mauroy

close all;
clear all;

% simulation parameters
t_end = 50;

x_interv = linspace(-1,2,100);
y_interv = linspace(-1,1,100);

% model parameters
a = 1; 
I = 0.05;
epsilon = 0.08;
gamma = 1;

param.a = a;
param.I = I;
param.epsilon = epsilon;
param.gamma = gamma;

f_dyn     =   @(t,X) FN(t,X,param);
options = odeset('RelTol',1e-12,'AbsTol',1e-300); 

pt_fix = roots([-1 (a+gamma) -(gamma*a+1) I]);
pt_fix(find(imag(pt_fix)~=0)) = [];
pt_fix = [pt_fix pt_fix];

J = [-3*pt_fix(1)^2+2*(a+gamma)*pt_fix(1)-a*gamma -1;epsilon -epsilon*gamma];
[u,v] = eig(J);
lambda1 = max([v(1,1) v(2,2)]);

%% computation of the Laplace averages

[x0,y0] = meshgrid(x_interv,y_interv);
point = [x0(:) y0(:)];

number_traject_per_loop = 100;
number_traject = length(point);
number_loops = ceil(number_traject/number_traject_per_loop);
average = zeros(1,numel(x0));

for test = 1 : number_loops;
    
    % setting the initial conditions
    index_low = 1+(test-1)*number_traject_per_loop;
    index_high = min(test*number_traject_per_loop,number_traject);
    number_traj_loop = index_high-index_low+1;
    init_cond = [point(index_low:index_high,1);point(index_low:index_high,2)];

    % ode integration
    [t,x] = ode45(f_dyn,[0 t_end],init_cond,options);
    
    % observable
    x1 = x(:,1:number_traj_loop);
    x2 = x(:,number_traj_loop+1:2*number_traj_loop);
    f = x1-pt_fix(1)+x2-pt_fix(2);
    
    % averaging
    average(index_low:index_high) = f(end,:).*(exp(-t(end)*lambda1)*ones(1,number_traj_loop));

end

data = reshape(average,size(x0));

%% plots

figure(1)
title('Isotables')
contour(x0,y0,abs(data),[0:0.2:5])

