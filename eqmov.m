% eqmov.m - Trabalho de Simulação - Parte 1
% Simulação numéricada dinâmica de voo logitudinal, não linear, de uma aeronave
% Alberto Romanhol Moreira - 2017051564
function [xderiv] = eqmov(x,u,param) 

%% DADOS FIXOS
m = 1200; % kg
g = 9.80665;
Iy = 1825; % kg m^2
S = 16; % m^2
c = 1.5; % m
rho = 1.05; % kg/m^3
CD0 = 0.027;
CDa = 1.121;
CLa_dot = 1.7;
Cma_dot = -7.27;
eT = 0; % rad

%% EQUAÇÕES
% Elementos do vetor de estado x
% x(1,i) Magnitude da velocidade m/s
% x(2,i) Ângulo de ataque rad
% x(3,i) Velocidade de arfagem rad/s
% x(4,i) Ângulo de arfagem rad

% Elemento do vetor de entrada u
% u(1,i) Deflexão de profundor rad
% u(2,i) Força propulsiva N

% Elementos do vetor de parâmetros ajustáveis param
% param(1) CL_0
% param(2) CL_alpha
% param(3) CL_q
% param(4) CL_deflexao_profundor
% param(5) Cm_0
% param(6) Cm_alpha
% param(7) Cm_q
% param(8) Cm_deflexao_profundor

[lin,tam] = size(x); % Tamanho do vetor x para o for
xderiv = zeros(lin,tam); % Definindo o tamnho do vetor de saida
for i = 1:tam
q_adm = (x(3,i)*c)/(2*x(1,i)); % q admensional
q_dash = 0.5*(x(1,i)^2)*rho*S; % Pressão dinâmica * Área

% CALCULO DE V PONTO
D = q_dash*(CD0 + CDa*x(2,i)); % Calculo de D
xderiv(1,i) = (1/m)*(-D - m*g*sin(x(4,i)-x(2,i)) + u(2,i)*cos(eT-x(2,i))); % V ponto

% CALCULO ALPHA PONTO
cl_term = (q_dash/(m*x(1,i)))*(param(1) + param(2)*x(2,i) + param(3)*q_adm + param(4)*u(1,i)); % Parâmetro que envolve os CL
qw_term = (1/(m*x(1,i)))*(-m*g*cos(x(4,i)-x(2,i)) - u(2,i)*sin(eT-x(2,i))); % Parâmetro que envolve qw
xderiv(2,i) = (x(3,i) - cl_term - qw_term)/((q_dash*CLa_dot*c)/(2*m*(x(1,i))^2) + 1); % Calculo de alpha ponto

a_dot_adm = (xderiv(2,i)*c)/(2*x(1,i)); % Alpha ponto admensional

% CALCULO Q PONTO
Cm = (param(5) + param(6)*x(2,i) + Cma_dot*a_dot_adm + param(7)*q_adm + param(8)*u(1,i)); % Calculo de Cm
M = q_dash*c*Cm; % Calculo do Momento
xderiv(3,i) = M/Iy; % q ponto

% TETA PONTO
xderiv(4,i) = x(3,i); % Teta ponto = q

end
end