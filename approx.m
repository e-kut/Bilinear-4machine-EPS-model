% Скрипт вычисляет матрицы A, N_k, B уравнения (23). Все формулы, на
% которые есть ссылки в этом коде, можно посмотреть в работе "Построение и
% анализ квадратичной аппроксимации нелинейной модели двухрайонной тестовой
% электроэнергетической системы в пространстве состояний".


% ШАГ 1: ВЫЧИСЛЕНИЕ МАТРИЦ ИЗ УРАВНЕНИЯ (20) И МАТРИЦЫ Y ПРОВОДИМОСТИ
clear all
sys_mat.init = fun_diagmat; % Матрицы Ad,Bd,Cd,Dd,Ad_2,Bd_2,Cd_2,Sd,Hd,B1,Y
[sys_mat.init.E,sys_mat.init.M,sys_mat.init.P] = fun_emp; % Матрицы E, M, P

% ШАГ 2: ВЫЧИСЛЕНИЕ МАТРИЦ ИЗ УРАВНЕНИЯ (22)
sys_mat.init = fun_sssys (sys_mat.init);

% ШАГ 3: ВЫЧИСЛЕНИЕ МАТРИЦЫ БИЛИНЕАРИЗОВАННОЙ СИСТЕМЫ (УРАВНЕНИЕ (23))
sys_mat.init = fun_bilinsys(sys_mat.init);

% ШАГ 4: УСТРАНЕНИЕ ЛИНЕЙНОЙ ЗАВИСИМОСТИ
sys_mat.lindep_sys = fun_lindep (sys_mat.init);

% ШАГ 5: УСТРАНЕНИЕ ИЗБЫТОЧНОСТИ В ВЕКТОРЕ ПЕРЕМЕННЫХ СОСТОЯНИЯ
sys_mat.reduce_sys = fun_reduce(sys_mat);

%% ОПЦИОНАЛЬНО: ПОСТРОЕНИЕ ГРАФИКОВ ПЕРЕХОДНЫХ ПРОЦЕССОВ
run ('trp_simu.m')

%% ============================ ФУНКЦИИ ================================

function [E,M,P] = fun_emp()
    Ei = struct;
    I = eye(6,6);
    for i = 1:4
        name_field = string(['Ei' num2str(i)]);
        tmpE = zeros(36,144);
        for k = 1:6
            l = 4 * (k-1) + i;
            tmpE(6*k-5:6*k,6*l-5:6*l) = I; 
        end
        Ei = setfield(Ei,name_field,tmpE);
    end
    E = blkdiag(Ei.Ei1,Ei.Ei2,Ei.Ei3,Ei.Ei4);

    % Матрица M
    Mi = struct;
    I = eye(2,2);
    for i = 1:4
        name_field = string(['Mi' num2str(i)]);
        tmpM = zeros(4,24);
        for k = 1:2
            l = 6 * (k-1) + i;
            tmpM(2*k-1:2*k,2*l-1:2*l) = I; 
        end
        Mi = setfield(Mi,name_field,tmpM);
    end
    M = [blkdiag(Mi.Mi1,Mi.Mi2,Mi.Mi3,Mi.Mi4) zeros(16,48);
        zeros(8,96) zeros(8,48)];

    % Матрица P
    Pi = struct;
    I = eye(6,6);
    for i = 1:4
        name_field = string(['Pi' num2str(i)]);
        tmpP = zeros(12,48);
        for k = 1:2
            l = 4 * (k-1) + i;
            tmpP(6*k-5:6*k,6*l-5:6*l) = I; 
        end
        Pi = setfield(Pi,name_field,tmpP);
    end
    P = [blkdiag(Pi.Pi1,Pi.Pi2,Pi.Pi3,Pi.Pi4) zeros(48,96)];
end

function sys_mat = fun_sssys(sys_mat)
    Y = sys_mat.Y;
    Ad = sys_mat.Ad;
    Ad_2 = sys_mat.Ad_2;
    Bd = sys_mat.Bd;
    Bd_2 = sys_mat.Bd_2;
    Cd = sys_mat.Cd;
    Cd_2 = sys_mat.Cd_2;
    Dd = sys_mat.Dd;
    Hd = sys_mat.Hd;
    Sd = sys_mat.Sd;
    E = sys_mat.E;
    M = sys_mat.M;
    P = sys_mat.P;
    B1 = sys_mat.B1;
    % Формируем символьный вектор x и x^(2) (x_kr)
    x_G1 = sym('x_g1', [6 1]);
    x_G2 = sym('x_g2', [6 1]);
    x_G3 = sym('x_g3', [6 1]);
    x_G4 = sym('x_g4', [6 1]);
    x = [x_G1; x_G2; x_G3; x_G4];
    x_kr = kron(x,x);
    % Единичная матрица
    I = eye(12);
    % Выражение для V в символьном виде (см. формулу (21))
    V = inv(Y - Dd - Sd * P * kron(I,x)) * (Cd * x + Cd_2 * E * kron(x,x));
    % Выражения для 2-го, 4-го и 5-го слагаемого формулы (20) в символьном виде
    term_2 = Bd * V;
    term_4 = Bd_2 * M * kron(V,V);
    term_5 = Hd * P * kron(V,x);
    % Раскладываем term_2,term_4,term_5 в ряд Тейлора до 2го порядка в точке x0 = 0
    TaylorSeries.term2.taylor = taylor(term_2,x,'Order',3);
    TaylorSeries.term4.taylor = taylor(term_4,x,'Order',3);
    TaylorSeries.term5.taylor = taylor(term_5,x,'Order',3);

    % Вычисление матриц из уравнения (22)

    % Матрицы B_D*B_lin и B_D*B_sqr
    sx = string(x);
    sx_kr = string(x_kr);
    tmp_tayseries = TaylorSeries.term2.taylor;
    Bd_Blin = zeros(size(tmp_tayseries,1),size(x,1));
    Bd_Bsqr = zeros(size(tmp_tayseries,1),size(x_kr,1));
    for i = 1:size(tmp_tayseries,1)
        [c, t] = coeffs(tmp_tayseries(i), x);
        st = string(t);
        if ~isempty(c)
            Bd_Blin (i,:) = get_lin(c,st,sx);
            Bd_Bsqr (i,:) = get_sqr(c,st,sx_kr);
        end 
    end
    % Матрицы B_D2*M*B_2lin и B_D2*M*B_2sqr
    tmp_tayseries = TaylorSeries.term4.taylor;
    Bd2M_B2lin = zeros(size(tmp_tayseries,1),size(x,1));
    Bd2M_B2sqr = zeros(size(tmp_tayseries,1),size(x_kr,1));
    for i = 1:size(tmp_tayseries,1)
        [c, t] = coeffs(tmp_tayseries(i), x);
        st = string(t);
        if ~isempty(c)
            Bd2M_B2lin (i,:) = get_lin(c,st,sx);
            Bd2M_B2sqr (i,:) = get_sqr(c,st,sx_kr);
        end 
    end
    % Матрицы H_D*P*H_lin и H_D*P*H_sqr
    tmp_tayseries = TaylorSeries.term5.taylor;
    HdP_Hlin = zeros(size(tmp_tayseries,1),size(x,1));
    HdP_Hsqr = zeros(size(tmp_tayseries,1),size(x_kr,1));
    for i = 1:size(tmp_tayseries,1)
        [c, t] = coeffs(tmp_tayseries(i), x);
        st = string(t);
        if ~isempty(c)
            HdP_Hlin (i,:) = get_lin(c,st,sx);
            HdP_Hsqr (i,:) = get_sqr(c,st,sx_kr);
        end 
    end
    % Матрицы А1, А2
    sys_mat.A1 = Ad + Bd_Blin + Bd2M_B2lin + HdP_Hlin;
    sys_mat.A2 = Ad_2 * E + Bd_Bsqr + Bd2M_B2sqr + HdP_Hsqr;
    % Перечень состояний
    sys_mat.x = string(x);
    sys_mat.x_kr = string(x_kr);
    sys_mat.x_full = [sys_mat.x;sys_mat.x_kr];
end

function c_lin = get_lin(c,t,x)
    % Достаём линейные коэффициенты
    c_lin = zeros (1,size(x,1));
    for i_x = 1 : size(x,1)
        id = find(strcmp(x(i_x),t));
        if ~isempty(id)
            c_lin(i_x) = c(id);
        end    
    end
end

function c_sqr = get_sqr (c,t,x_kr)
    % Достаём квадратичные коэффициенты
    exception_x_kr = string;
    c_sqr = zeros (1,size(x_kr,1));
    for i_x = 1 : size(x_kr,1)
        if sum(x_kr(i_x) == exception_x_kr) == 1 % если элемент из x_kr уже был
            continue %  то переходим к следующему шагу
        end
        tmp_x = x_kr(i_x);
        tf = strcmp(tmp_x,t);
        id = find(tf);
        if ~isempty(id)
            c_sqr(i_x) = c(id);
        end
        exception_x_kr(i_x) = tmp_x;
    end
end

function sys_mat = fun_bilinsys(sys_mat)
    A1 = sys_mat.A1;
    A2 = sys_mat.A2;
    B1 = sys_mat.B1;
    % Матрица А билинеаризованной системы
    n = size(A1,1);
    n2 = n + n * n;
    I = eye(n, n);
    A21 = kron(A1, I) + kron(I, A1);
    A = zeros(n + n^2, n + n^2);
    A(1:n, 1:n) = A1;
    A(1:n, n+1:n2) = A2;
    A(n+1:n2, n+1:n2) = A21;
    sys_mat.A21 = A21;
    sys_mat.A = A;
    % Матрица B (управления) билинеаризованной системы
    B = zeros(size(A,1), 8);
    B(1:n,:) = B1;
    sys_mat.B = B;
    % Матрица N билинеаризованной системы
    for i = 1: size(B1,2)
        eval(['B20',num2str(i), ' = kron(B1(:,i), I) + kron(I, B1(:,i));'])
        eval(['sys_mat.N.N',num2str(i), ' = zeros(size(A));'])
        eval(['sys_mat.N.N',num2str(i),'(n+1:n2, 1:n)' ' = B20',num2str(i),';'])
    end
end

function lindep_sys = fun_lindep (sys_mat)
    A1 = sys_mat.A1; % линейная часть матрицы билинейной системы
    A2 = sys_mat.A2; % квадратичная часть матрицы билинейной системы
    B1 = sys_mat.B1;
    N = sys_mat.N;
    % вычитаем строки:
    A1([7,8],:) = A1([7,8],:) - A1([1,2],:); 
    A1([13,14],:) = A1([13,14],:) - A1([1,2],:);
    A1([19,20],:) = A1([19,20],:) - A1([1,2],:);
    A1([1,2],:) = A1([1,2],:) - A1([1,2],:);
    A2([7,8],:) = A2([7,8],:) - A2([1,2],:); 
    A2([13,14],:) = A2([13,14],:) - A2([1,2],:);
    A2([19,20],:) = A2([19,20],:) - A2([1,2],:);
    A2([1,2],:) = A2([1,2],:) - A2([1,2],:);
    % перестраиваем матрицу A21 и матрицу А билинеаризованной системы
    n = size(A1,1);
    n2 = n + n * n;
    I = eye(n, n);
    A21 = kron(A1, I) + kron(I, A1);
    A_lnd = zeros(n + n^2, n + n^2);
    A_lnd(1:n, 1:n) = A1;
    A_lnd(1:n, n+1:n2) = A2;
    A_lnd(n+1:n2, n+1:n2) = A21;
    % матрицы линейной и квадратичной аппроксимации
    lindep_sys.A1 = A1;
    lindep_sys.A = A_lnd;
    % матрицы B и N:
    % вычитаем строки:
    B1([7,8],:) = B1([7,8],:) - B1([1,2],:); 
    B1([13,14],:) = B1([13,14],:) - B1([1,2],:);
    B1([19,20],:) = B1([19,20],:) - B1([1,2],:);
    B1([1,2],:) = B1([1,2],:) - B1([1,2],:);
    % матрица B (управления) билинеаризованной системы
    B = zeros(size(A_lnd,1), 8);
    B(1:n,:) = B1;
    % Матрица N билинеаризованной системы
    for i = 1: size(B1,2)
        eval(['B20',num2str(i), ' = kron(B1(:,i), I) + kron(I, B1(:,i));'])
        eval(['N.N',num2str(i), ' = zeros(size(A_lnd));'])
        eval(['N.N',num2str(i),'(n+1:n2, 1:n)' ' = B20',num2str(i),';'])
    end
    lindep_sys.B1 = B1;
    lindep_sys.B = B;
    lindep_sys.N = N;
end

function reduce_sys = fun_reduce(sys_mat)
    A = sys_mat.lindep_sys.A;
    B = sys_mat.lindep_sys.B;
    X = sys_mat.init.x_full;
    s = 0;
    r = 0;
    ID = {};
    for p = 1:size(A,1)
        tmp = find(X(p) == X);
        ID{p,1} = tmp;
        if tmp(1) == p
            s = s + 1;
            id(s,1) = p; % список номеров переменных без повторений
        else
            r = r + 1;
            rmv(r,1) = p; % список номеров повторяющихся переменных, подлежащих удалению
            svd(r,1) = tmp(1); % список номеров повторяющихся переменных, подлежащих сохранению
        end
    end
    Ared = A;
    Ared(:,svd) = A(:,svd) + A(:,rmv);
    Ared = Ared(id,:);
    Ared = Ared(:,id);
    Bred = B(id,:);
    Xred = X(id);
    N = struct;
    for k = 1:size(B,2)
        name_field = string(['N' num2str(k)]);
        tN = getfield(sys_mat.lindep_sys.N,name_field);
        tN(:,svd) = tN(:,svd) + tN(:,rmv);
        Nred = tN(id,:);
        Nred = Nred(:,id);
        N = setfield(N,name_field,Nred);
    end
    reduce_sys.A = Ared;
    reduce_sys.A1 = Ared(1:24,1:24);
    reduce_sys.B = Bred;
    reduce_sys.B1 = Bred(1:24,:);
    reduce_sys.N = N;
    reduce_sys.x_full = Xred;
    reduce_sys.ids = id;
end