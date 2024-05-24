function FN=FN(t,x,param)

a=param.a;
I=param.I;
epsilon=param.epsilon;
gamma=param.gamma;


if size( x, 2) >1
    x=x(:);
end

nPoints = length( x) / 2;
ix = (0 * nPoints + 1) : (1 * nPoints);
iy = (1 * nPoints + 1) : (2 * nPoints);

FN=[-x(iy)-x(ix).*(x(ix)-1).*(x(ix)-a)+I;epsilon*(x(ix)-gamma*x(iy))];
