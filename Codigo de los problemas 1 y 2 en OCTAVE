clc; clear;

%  PROBLEMA 1: GAUSS CON PIVOTEO Y REDONDEO A 3 DÍGITOS

disp("--- Problema 1: Gauss con Redondeo a 3 dígitos ---");
format short g; % Formato para ver 3 dígitos significativos

% --- Función Auxiliar de Redondeo a 'sigs' dígitos significativos ---
function r = round_sig(val, sigs)
    if val == 0
        r = 0;
    else
        % Formatea a 'sigs' dígitos (ej. '%.3g') y lo vuelve a leer
        r = sscanf(sprintf(['%.' num2str(sigs) 'g'], val), '%f');
    end
endfunction

% --- Configuración del Problema 1 ---
n_sig = 3; % Dígitos significativos [fuente: 6]
% Matriz aumentada [A|b] (con corrección en Ec. 3 [fuente: 7])
Ab = [1.19,  2.11, -100,   1.22;
      14.2, -0.122, 12.2,   3.44;
      0,     100,   -99.9,  2.15];
n = 3;

% --- Fase de Eliminación (Silenciosa) ---
for k = 1:n-1
    % 1. Pivoteo Parcial
    [max_val, max_idx] = max(abs(Ab(k:n, k)));
    max_idx = max_idx + (k-1);
    if max_idx != k
        temp_row = Ab(k, :);
        Ab(k, :) = Ab(max_idx, :);
        Ab(max_idx, :) = temp_row;
    end
    % 2. Eliminación con redondeo
    for i = k+1:n
        multiplicador = round_sig(Ab(i, k) / Ab(k, k), n_sig);
        for j = k:n+1
            term = round_sig(multiplicador * Ab(k, j), n_sig);
            Ab(i, j) = round_sig(Ab(i, j) - term, n_sig);
        end
    end
end

% --- Fase de Sustitución Hacia Atrás (Silenciosa) ---
x1 = zeros(n, 1);
for i = n:-1:1
    sum = 0;
    for j = i+1:n
        term = round_sig(Ab(i, j) * x1(j), n_sig);
        sum = round_sig(sum + term, n_sig);
    end
    numerador = round_sig(Ab(i, n+1) - sum, n_sig);
    x1(i) = round_sig(numerador / Ab(i, i), n_sig);
end

disp("Solución (Gauss con redondeo):");
fprintf("x1 = %f\nx2 = %f\nx3 = %f\n", x1(1), x1(2), x1(3));


%% =================================================================
%  PROBLEMA 2: MÉTODOS ITERATIVOS (ALTA PRECISIÓN)
%  =================================================================
disp(" "); % Salto de línea
disp("--- Problema 2: Métodos Iterativos (Alta Precisión) ---");
format long; % Formato de alta precisión

% --- Configuración del Problema 2 --- [fuente: 9, 10, 11, 12]
A = [5, -1,  2;
    -1,  4,  1;
     2,  1,  5];
b = [12; -2; 10];

tol = 1e-3; % Tolerancia [fuente: 15]
x_inicial = [0; 0; 0];
max_iter = 100;

% Descomposición A = D + L + U
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);

% --- Análisis de Convergencia ---
disp(" ");
disp("Análisis de Convergencia:");
T_J = -inv(D) * (L + U); % Matriz Jacobi [fuente: 13]
rho_J = max(abs(eig(T_J))); % Radio espectral Jacobi [fuente: 14]
T_G = -inv(D + L) * U; % Matriz Gauss-Seidel [fuente: 16]
rho_G = max(abs(eig(T_G))); % Radio espectral G-S [fuente: 16]

fprintf("Radio Espectral rho(Jacobi) = %f\n", rho_J);
fprintf("Radio Espectral rho(Gauss-Seidel) = %f\n", rho_G);
fprintf("rho(G-S) < rho(J), por lo tanto G-S converge más rápido.\n\n");

% --- c. Iteraciones de Jacobi --- [fuente: 15]
disp("Iteraciones de Jacobi (k, x1, x2, x3, Error):");
x = x_inicial;
error = inf;
k = 0;
while error > tol && k < max_iter
    x_anterior = x;
    x(1) = (b(1) - A(1,2)*x_anterior(2) - A(1,3)*x_anterior(3)) / A(1,1);
    x(2) = (b(2) - A(2,1)*x_anterior(1) - A(2,3)*x_anterior(3)) / A(2,2);
    x(3) = (b(3) - A(3,1)*x_anterior(1) - A(3,2)*x_anterior(2)) / A(3,3);
    error = norm(x - x_anterior, inf);
    k = k + 1;
    fprintf("%2d | %11.7f | %11.7f | %11.7f | %11.7f\n", k, x(1), x(2), x(3), error);
end
fprintf("Solución (Jacobi) en %d iteraciones:\n", k);
disp(x);

% --- d. Iteraciones de Gauss-Seidel --- [fuente: 16]
disp(" ");
disp("Iteraciones de Gauss-Seidel (k, x1, x2, x3, Error):");
x = x_inicial;
error = inf;
k = 0;
while error > tol && k < max_iter
    x_anterior = x;
    x(1) = (b(1) - A(1,2)*x_anterior(2) - A(1,3)*x_anterior(3)) / A(1,1);
    x(2) = (b(2) - A(2,1)*x(1)          - A(2,3)*x_anterior(3)) / A(2,2);
    x(3) = (b(3) - A(3,1)*x(1)          - A(3,2)*x(2))           / A(3,3);
    error = norm(x - x_anterior, inf);
    k = k + 1;
    fprintf("%2d | %11.7f | %11.7f | %11.7f | %11.7f\n", k, x(1), x(2), x(3), error);
end
fprintf("Solución (Gauss-Seidel) en %d iteraciones:\n", k);
disp(x);
