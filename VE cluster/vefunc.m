function SS=vefunc(y)
global k1 k2 k3 k4 k7 k9 k10 c0 cmin

SS =[-k1*y(1)*y(2)+k2*y(2)-2*k3*y(1)^2+k4*y(1)*(c0-y(3))/(c0-y(3)+cmin);
     -3*k1*y(1)*y(2)-2*k2*y(2)-k3*y(1)^2+k7*y(4)+k9*y(3);
     2*k4*y(1)*(c0-y(3))/(c0-y(3)+cmin)-k9*y(3)-k10*y(3);
     2*k1*y(1)*y(2)+k2*y(2)+k3*y(1)^2-k7*y(4)];