% mech 6481 aeroelasticity pk method
clc;
clear;
close all;

%% code evaluation 
% variations 
a = -1/5;
e = -1/10;
xtheta = e-a;
mu = 20;
rs = 6/25; 
sigma = 2/5;

N = 101;


V_vec = linspace(0,4,N);
r1_v = zeros(1,N);
r2_v = zeros(1,N); 
r3_v = zeros(1,N);
r4_v = zeros(1,N); 
K_v = zeros(4,N);
Gamma = zeros(4,N);
omega_Dimensionless = zeros(4,N);
k_guess = 0.1;
Tolerance_v = 1e-7;

  for s=2:N
      V = V_vec (s);
      k = k_guess;

      for q = 1 : 4
         k_distance_v = 1;
          while abs(k_distance_v)>Tolerance_v
              C_fun_v = (0.01365+(0.2808*1i*k)-(k^2/2))/(0.01365+(0.3455*1i*k)-(k^2));
              f11_v = (sigma/V)^2-(k^2/mu)+((2*1i*k*C_fun_v)/mu);
              f12_v = ((k*(1i+(a*k)))+(C_fun_v*(2+(1i*k)*(1-(2*a)))))/mu;
              f21_v = ((a*k^2)-(1i*k*C_fun_v*(1+(2*a))))/mu;
              f22_v =((8*mu*rs/V^2)+(4*1i*(1+2*a))*((2*1i-k*(1-2*a)))*C_fun_v-k*(k-4*1i+8*a*(1i+a*k)))/(8*mu);
              poly_v = [(rs-xtheta^2) (0) (f22_v+(f11_v*rs)-(f21_v*xtheta)-(f12_v*xtheta)) (0) ((f11_v*f22_v)-(f12_v*f21_v))];
              rts_v = roots(poly_v);
              [not_v,arrange_v] = sort(imag(rts_v));
              rts_1_v = rts_v(arrange_v);
              k_new_v = imag(rts_1_v(q));
              Gamma_new_v = real(rts_1_v(q))/imag(rts_1_v(q));
              k_distance_v = k_new_v - k;
              k = k_new_v;
          end
          K_v(q,s) = k_new_v;
          Gamma(q,s) = Gamma_new_v;
          omega_Dimensionless(q,s) = k_new_v*V;
      end
  end

figure (1);
plot(V_vec,omega_Dimensionless(1,:),'b*','MarkerSize',2);
hold on;
plot(V_vec,omega_Dimensionless(2,:),'ko','MarkerSize',2);
plot(V_vec,omega_Dimensionless(3,:),'rd','MarkerSize',2);
plot(V_vec,omega_Dimensionless(4,:),'gs','MarkerSize',2);
grid on;
xlabel('V= U/(b \omega_\theta)');
ylabel('(\omega/\omega_\theta)');
figure (2);
plot(V_vec,Gamma(1,:),'b*','MarkerSize',2);
hold on;
plot(V_vec,Gamma(2,:),'ko','MarkerSize',2);
plot(V_vec,Gamma(3,:),'rd','MarkerSize',2);
plot(V_vec,Gamma(4,:),'gs','MarkerSize',2);
grid on;
xlabel('V= U/(b \omega_\theta)');
ylabel('\gamma');