ms = 1250;
ksl = 50000;
ksr = 50000;
bsl = 3000;
bsr = 3000;
mu = 150;
ktl = 350000;
ktr = 350000;
t = 2.5;
h1 = 1;
h2 = 2.2;
Is = 2000;
Iu = 200;
g = 9.81;

rho = [16, 16, 40000, 16, 10^-9, 10^-9];

Iu_ = Iu + mu*h1^2;
Is_ = Is + ms*h2^2;

A_1 = [0, 1, 0, 0];
A_2 = [-(ktl+ktr+ksl+ksr)*t^2/4 + mu*g*h1, -(bsl+bsr)*t^2/4, (ksl+ksr)*t^2/4, (bsl+bsr)*t^2/4] / Iu_;
A_3 = [0, 0, 0, 1];
A_4 = [(ksl+ksr)*t^2/4, (bsl+bsr)*t^2/4, -(ksl+ksr)*t^2/4 + ms*g*h2, -(bsl+bsr)*t^2/4] / Is_;
A = [A_1; A_2; A_3; A_4];

B = [0, 0; -t/(2*Iu_), t/(2*Iu_); 0, 0; t/(2*Is_), -t/(2*Is_)];

L = [0; mu*h1/Iu_; 0; ms*h2/Is_];

Q = diag(rho(1:4));
R = diag(rho(5:6));
N = 0;


[K, S, e] = lqr(A, B, Q, R, N);


%% System Response

ay0 = 1;
t = 0:0.01:25;

in_pulse = zeros(1, length(t));
in_pulse(abs(t-1.2)<0.2) = ay0*(1 - 1/0.2*abs(t(abs(t-1.2)<0.2)-1.2));

in_rect = zeros(1, length(t));
in_rect(find(abs(t-1)<=0.001):end) = ay0;
in_rect(find(abs(t-2)<=0.001):end) = 0;


[s_a, s_b] = ss2tf(A, L, [1, 0, 0, 0], 0);
ta1_p = tf(s_a, s_b);

[s_a, s_b] = ss2tf(A, L, [0, 0, 1, 0], 0);
ta2_p = tf(s_a, s_b);

[s_a, s_b] = ss2tf(A-B*K, L, [1, 0, 0, 0], 0);
ta1_a = tf(s_a, s_b);

[s_a, s_b] = ss2tf(A-B*K, L, [0, 0, 1, 0], 0);
ta2_a = tf(s_a, s_b);


figure(1);
y = lsim(ta1_p, in_pulse, t);
plot(t, y);
hold on;
y = lsim(ta1_a, in_pulse, t);
plot(t, y)
title('Roll Angle of Unsprung Mass (Impulse Input)');
xlabel('time (s)');
ylabel('\alpha_1 (rad)');
legend({'Passive', 'Active'});

figure(2);
y = lsim(ta2_p, in_pulse, t);
plot(t, y);
hold on;
y = lsim(ta2_a, in_pulse, t);
title('Roll Angle of Sprung Mass (Impulse Input)');
plot(t, y);
xlabel('time (s)');
ylabel('\alpha_2 (rad)');
legend({'Passive', 'Active'});


figure(3);
y = lsim(ta1_p, in_rect, t);
plot(t, y);
hold on;
y = lsim(ta1_a, in_rect, t);
plot(t, y)
title('Roll Angle of Unsprung Mass (Rectangular Input)');
xlabel('time (s)');
ylabel('\alpha_1 (rad)');
legend({'Passive', 'Active'});

figure(4);
y = lsim(ta2_p, in_rect, t);
plot(t, y);
hold on;
y = lsim(ta2_a, in_rect, t);
title('Roll Angle of Sprung Mass (Rectangular Input)');
plot(t, y);
xlabel('time (s)');
ylabel('\alpha_2 (rad)');
legend({'Passive', 'Active'});