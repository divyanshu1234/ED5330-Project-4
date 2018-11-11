ms = 300;
mu = 50;
ks = 15000;
kt = 150000;
bs = 900;
wn = {6.5*10^-2, 57.5*10*2};

s = tf('s');

d = mu*ms*s^4 + (mu + ms)*bs*s^3 + (mu*ks + ms*ks + ms*kt)*s^2 + bs*kt*s + ks*kt;
Ta = (kt*s * (bs*s + ks)) / d;
Tr = (-kt*ms*s) / d;
Tt = -(mu*ms*s^3 + (mu + ms)*bs*s^2 + (mu + ms)*ks * s) / d;

figure(1);
bode(Ta, wn);
title('Bode - Acceleration Transfer Function');
% margin(Ta);

figure(2);
bode(Tr, wn);
title('Bode - Rattle Space Transfer Function');
% margin(Tr);

figure(3);
bode(Tt, wn);
title('Bode - Tyre Deflection Transfer Function');
% margin(Tt);