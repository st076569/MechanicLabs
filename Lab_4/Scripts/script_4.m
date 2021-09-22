% Крутильные колебания вала с дисками (Баталов Семен, 351)
inputData = readmatrix('input_data.csv');

% Радиусы дисков
r_1 = inputData(1, 2);
r_2 = inputData(1, 3);
r_3 = inputData(1, 4);

% Толщины дисков
d_1 = inputData(1, 5);
d_2 = inputData(1, 6);
d_3 = inputData(1, 7);

% Расстояния между дисками
l_1 = inputData(1, 8);
l_2 = inputData(1, 9);

% Радиус вала и длина эксцентрика
r = inputData(1, 10);
e = inputData(1, 11);

% Известные константы
rho = inputData(1, 12);
g   = inputData(1, 13);
c_p = inputData(1, 14);

% Моменты инерции дисков
I_1 = momInertion(rho, r_1, d_1);
I_2 = momInertion(rho, r_2, d_2);
I_3 = momInertion(rho, r_3, d_3);

% Жесткость вала на скручивание
c_1 = torsionStiffness(g, r, l_1);
c_2 = torsionStiffness(g, r, l_2);
c_3 = c_1 + c_2 + 2 * c_p * r_2.^2;

% Коэффициенты многочлена
a_0 = I_1 * I_2 * I_3;
a_1 = - c_2 * I_1 * I_2 - c_3 * I_1 * I_3 - c_1 * I_2 * I_3;
a_2 = c_2 * (c_3 - c_2) * I_1 + c_1 * c_2 * I_2 + c_1 * (c_3 - c_1) * I_3;
a_3 = - c_1 * c_2 * (c_3 - c_1 - c_2);

% Решаем уравнение
syms y;
eqn = a_0 * y^3 + a_1 * y^2 + a_2 * y + a_3 == 0;
s = solve(eqn, y);

% Искомые частоты собственных колебаний системы
omega_1 = sqrt(s(1));
omega_2 = sqrt(s(2));
omega_3 = sqrt(s(3));

% Создаем матрицы системы для omega_1, omega_2, omega_3 и столбец 
% свободных членов
b = [0; r_2 * c_p * e; 0];
M_1 = createMatrix(c_1, c_2, c_3, I_1, I_2, I_3, omega_1);
M_2 = createMatrix(c_1, c_2, c_3, I_1, I_2, I_3, omega_2);
M_3 = createMatrix(c_1, c_2, c_3, I_1, I_2, I_3, omega_3);

% Находим главные формы собственных колебаний и нормируем их
phi_omega_1 = [detNum(M_1, b, 1); detNum(M_1, b, 2); detNum(M_1, b, 3)];
phi_omega_2 = [detNum(M_2, b, 1); detNum(M_2, b, 2); detNum(M_2, b, 3)];
phi_omega_3 = [detNum(M_3, b, 1); detNum(M_3, b, 2); detNum(M_3, b, 3)];
phi_omega_1 = phi_omega_1 / norm(phi_omega_1);
phi_omega_2 = phi_omega_2 / norm(phi_omega_2);
phi_omega_3 = phi_omega_3 / norm(phi_omega_3);

% Вывод данных на экран
fprintf('\n > Входные данные : ');
fprintf('\n r_1 : %.5f', r_1);
fprintf('\n r_2 : %.5f', r_2);
fprintf('\n r_3 : %.5f', r_3);
fprintf('\n d_1 : %.5f', d_1);
fprintf('\n d_2 : %.5f', d_2);
fprintf('\n d_3 : %.5f', d_3);
fprintf('\n l_1 : %.5f', l_1);
fprintf('\n l_2 : %.5f', l_2);
fprintf('\n r   : %.5f', r);
fprintf('\n e   : %.5f', e);
fprintf('\n rho : %.5f', rho);
fprintf('\n g   : %.5f', g);
fprintf('\n c_p : %.5f', c_p);
fprintf('\n');
fprintf('\n > Моменты инерции дисков : ');
fprintf('\n I_1 : %.5f', I_1);
fprintf('\n I_2 : %.5f', I_2);
fprintf('\n I_3 : %.5f', I_3);
fprintf('\n');
fprintf('\n > Жескости на скручивание : ');
fprintf('\n c_1 : %.5f', c_1);
fprintf('\n c_2 : %.5f', c_2);
fprintf('\n c_3 : %.5f', c_3);
fprintf('\n');
fprintf('\n > Коэффициенты многочлена : ');
fprintf('\n a_0 : %.5f', a_0);
fprintf('\n a_1 : %.5f', a_1);
fprintf('\n a_2 : %.5f', a_2);
fprintf('\n a_3 : %.5f', a_3);
fprintf('\n');
fprintf('\n > Искомые угловые частоты : ');
fprintf('\n omega_1 : %.5f', omega_1);
fprintf('\n omega_2 : %.5f', omega_2);
fprintf('\n omega_3 : %.5f', omega_3);
fprintf('\n');
fprintf('\n > Искомые формы колебаний : ');
fprintf('\n phi_omega_1(1) : %.5f', phi_omega_1(1));
fprintf('\n phi_omega_1(2) : %.5f', phi_omega_1(2));
fprintf('\n phi_omega_1(3) : %.5f', phi_omega_1(3));
fprintf('\n');
fprintf('\n phi_omega_2(1) : %.5f', phi_omega_2(1));
fprintf('\n phi_omega_2(2) : %.5f', phi_omega_2(2));
fprintf('\n phi_omega_2(3) : %.5f', phi_omega_2(3));
fprintf('\n');
fprintf('\n phi_omega_3(1) : %.5f', phi_omega_3(1));
fprintf('\n phi_omega_3(2) : %.5f', phi_omega_3(2));
fprintf('\n phi_omega_3(3) : %.5f', phi_omega_3(3));
fprintf('\n');

% Момент инерции диска
function I = momInertion(rho, r, d)
    I = 0.5 * rho * pi * r.^4 * d;
end

% Жесткость вала на скручивание
function c = torsionStiffness(g, r, l)
    c = 0.5 * g * pi * r.^4 / l;
end

% Матрица системы
function m = createMatrix(c_1, c_2, c_3, I_1, I_2, I_3, omega)
    m = [c_1 - omega.^2 * I_1 -c_1 0;
         -c_1 c_3 - omega.^2 * I_2 -c_2;
         0 -c_2 c_2 - omega.^2 * I_3];
end

% Определитель матрицы системы с измененным столбцом
function f = detNum(m, b, num)
    new_m = m;
    new_m(1:3, num) = b;
    f = det(new_m);
end
