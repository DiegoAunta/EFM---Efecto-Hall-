% Definición de constantes
e = 1.602E-19;  % Carga elemental del electrón [C]
kb = 1.38E-23;  % Constante de Boltzmann [J/K]

% Definir incertidumbre
dT = 0.05;  % Incertidumbre de la temperatura
dVH = 0.01;  % Ejemplo de incertidumbre en VH, ajusta según tus datos

% Parámetros iniciales
parametros_iniciales = [1.17e21, 1.99e26, 0.74, 1.81];

% Función para el ajuste
VoltajeTemp = @(params, x) ((params(1) + (sqrt(params(1)^2 / 4 + params(2)^2 * exp(-params(3) * 13025.9 ./ x)) - params(1) / 2) * (1 - params(4)^2)) ./ ((params(1) + (sqrt(params(1)^2 / 4 + params(2)^2 * exp(-params(3) * 13025.9 ./ x)) - params(1) / 2) * (1 + params(4))) .^ 2) * 6.41e19);

% Función de ajuste para nlinfit
modelo = @(params, x) VoltajeTemp(params, x);

% Ajuste de curva
[parametros,~,~, covarianza] = nlinfit(T, VH, modelo, parametros_iniciales, 'Weights', 1 ./ (dVH * ones(size(VH))));

% Obtener parámetros y sus incertidumbres
A = parametros(1);
B = parametros(2);
C = parametros(3);
D = parametros(4);

incertidumbres = sqrt(diag(covarianza));
dA = incertidumbres(1);
dB = incertidumbres(2);
dC = incertidumbres(3);
dD = incertidumbres(4);

% Graficar resultados
figure;
hold on
plot(T, VoltajeTemp(parametros, T),LineWidth=1, Color="b");
% plot(T, VH,LineStyle="none",Marker=".",MarkerSize=20,Color="b");
xlabel('T[K]');
ylabel('V_H[V]');
set(gca, "Fontsize", 30, "FontName", "Cambria Math") % Opciones de fuente y tamaño
grid on
grid minor

% Cambiar la escala del eje Y a notación científica
ax = gca; % Obtén el objeto de ejes actual
% ax.XAxis.Exponent = -3;
ax.YAxis.Exponent = -3;

% Matriz para exportar
M = [A, dA;B, dB;C, dC;D, dD];

% Exportar la matriz a un archivo de texto
writematrix(M, 'VHvsT.txt', 'Delimiter', 'tab');