% Критические угловые скорости вращения вала (Баталов Семен, 351)
inputData = readmatrix('input_data.csv');

% Входные данные
d_0 = [inputData(1, 2) inputData(2, 2)];
l_0 = [inputData(1, 3) inputData(2, 3)];
l_1 = [inputData(1, 4) inputData(2, 4)];
m_1 = [inputData(1, 5) inputData(2, 5)];
E = [inputData(1, 6) inputData(2, 6)];

% Линейная плотность и момент инерции сечения вала
p = compute_p(m_1, l_1);
I = compute_I(d_0);

% Расчет критических угловых скоростей
omega_1 = compute_omega(1, l_0, p, I, E);
omega_2 = compute_omega(2, l_0, p, I, E);
omega_3 = compute_omega(3, l_0, p, I, E);

% Вывод данных на экран
fprintf('\n > Входные данные (значения) : ');
fprintf('\n d_0       : %.5f', d_0(1));
fprintf('\n l_0       : %.5f', l_0(1));
fprintf('\n l_1       : %.5f', l_1(1));
fprintf('\n m_1       : %.5f', m_1(1));
fprintf('\n E         : %.5f', E(1));
fprintf('\n > Входные данные (погрешности) : ');
fprintf('\n delta_d_0 : %.5f', d_0(2));
fprintf('\n delta_l_0 : %.5f', l_0(2));
fprintf('\n delta_l_1 : %.5f', l_1(2));
fprintf('\n delta_m_1 : %.5f', m_1(2));
fprintf('\n');
fprintf('\n > Вспомогательные величины (значения) : ');
fprintf('\n p       : %.5f', p(1));
fprintf('\n I       : %.12f', I(1));
fprintf('\n > Вспомогательные величины (погрешности) : ');
fprintf('\n delta_p : %.5f', p(2));
fprintf('\n delta_I : %.12f', I(2));
fprintf('\n');
fprintf('\n > Искомые критические частоты (значения) : ');
fprintf('\n omega_1        : %.5f', omega_1(1));
fprintf('\n omega_2        : %.5f', omega_2(1));
fprintf('\n omega_3        : %.5f', omega_3(1));
fprintf('\n > Искомые критические частоты (погрешности) : ');
fprintf('\n delta_omega_1  : %.5f', omega_1(2));
fprintf('\n delta_omega_2  : %.5f', omega_2(2));
fprintf('\n delta_omega_3  : %.5f', omega_3(2));
fprintf('\n');

% Линейная плотность вала
function f = compute_p(m_1, l_1)
    % Символьное выражение для p
    syms p(m, l)
    p = m / l;
    
    % Расчет частных производных
    der_m = double(subs(diff(p, m), {m, l}, [m_1(1) l_1(1)]));
    der_l = double(subs(diff(p, l), {m, l}, [m_1(1) l_1(1)]));
    
    % Косвенная погрешность измерения p
    delta_p = sqrt((der_m * m_1(2))^2 + (der_l * l_1(2))^2);
    f = [double(subs(p, {m, l}, [m_1(1) l_1(1)])) delta_p];
end

% Момент инерции площади поперечного сечения вала
function f = compute_I(d_1)
    % Символьное выражение для I
    syms I(d)
    I = pi * d^4 / 64;
    
    % Расчет частных производных
    der_I = double(subs(diff(I, d), {d}, [d_1(1)]));
    
    % Косвенная погрешность измерения I
    delta_I = sqrt((der_I * d_1(2))^2);
    f = [double(subs(I, {d}, [d_1(1)])) delta_I];
end

% Критическая угловая скорость
function f = compute_omega(n_1, l_1, p_1, I_1, E_1)
    % Символьное выражение для omega
    syms omega(n, l, p, I, E)
    omega = (pi^2 * n^2 / l^2) * sqrt(E * I / p);
    
    % Расчет частных производных
    der_l = double(subs(diff(omega, l), {n, l, p, I, E}, ...
        [n_1 l_1(1) p_1(1) I_1(1) E_1(1)]));
    der_p = double(subs(diff(omega, p), {n, l, p, I, E}, ...
        [n_1 l_1(1) p_1(1) I_1(1) E_1(1)]));
    der_I = double(subs(diff(omega, I), {n, l, p, I, E}, ...
        [n_1 l_1(1) p_1(1) I_1(1) E_1(1)]));
    
    % Косвенная погрешность измерения omega
    delta_omega = sqrt((der_l * l_1(2))^2 + (der_p * p_1(2))^2 + ...
        (der_I * I_1(2))^2);
    f = [double(subs(omega, {n, l, p, I, E}, ...
        [n_1 l_1(1) p_1(1) I_1(1) E_1(1)])) delta_omega];
end