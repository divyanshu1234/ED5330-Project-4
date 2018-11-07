ms = 300;
mu = 50;
ks = 15000;
kt_list = [100000, 150000, 200000];
bs = 900;
wn = {6.5*10^-2, 57.5*10*2};

gm_pm = [];

s = tf('s');

for kt = kt_list
    d = mu*ms*s^4 + (mu + ms)*bs*s^3 + (mu*ks + ms*ks + ms*kt)*s^2 + bs*kt*s + ks*kt;
    Ta = (kt*s * (bs*s + ks)) / d;
    Tr = (-kt*ms*s) / d;
    Tt = -(mu*ms*s^3 + (mu + ms)*bs*s^2 + (mu + ms)*ks * s) / d;

    figure(1);
    bode(Ta, wn);
    title('Bode - Acceleration Transfer Function');
    hold on;
    [gm, pm] = margin(Ta);
    gm_pm = [gm_pm; gm, pm];

    figure(2);
    bode(Tr, wn);
    title('Bode - Rattle Space Transfer Function');
    hold on;
    [gm, pm] = margin(Tr);
    gm_pm = [gm_pm; gm, pm];

    figure(3);
    bode(Tt, wn);
    title('Bode - Tyre Deflection Transfer Function');
    hold on;
    [gm, pm] = margin(Tt);
    gm_pm = [gm_pm; gm, pm];
end

figure(1); legend({'kt=100000 N/m', 'bs=150000 N/m', 'bs=200000 N/m'});
figure(2); legend({'kt=100000 N/m', 'bs=150000 N/m', 'bs=200000 N/m'});
figure(3); legend({'kt=100000 N/m', 'bs=150000 N/m', 'bs=200000 N/m'});