% This script calculate the matrices A, N_k, B from equation (23). All
% formulas referenced in this code can be found in the doc.pdf


% STEP#1: CALCULATION OF MATRICES FROM EQUATION (20) & ADDMITANCE MATRIX (Y)
clear all
sys_mat.init = fun_diagmat; % Matrices Ad,Bd,Cd,Dd,Ad_2,Bd_2,Cd_2,Sd,Hd,B1,Y
[sys_mat.init.E,sys_mat.init.M,sys_mat.init.P] = fun_emp; % Matrices E,M,P

% STEP#2: CALCULATION OF MATRICES FROM EQUATION (22)
sys_mat.init = fun_sssys (sys_mat.init);

% STEP#3: CALCULATION OF BILINEAR SYSTEM MATRICES (EQUATION (23))
sys_mat.init = fun_bilinsys(sys_mat.init);

% STEP#4: ELIMINATION OF LINEAR DEPENDENCE
sys_mat.lindep_sys = fun_lindep (sys_mat.init);

% STEP 5: ELIMINATING REDUNDANCY IN THE VECTOR OF STATE VARIABLES
sys_mat.reduce_sys = fun_reduce(sys_mat);

%% OPTIONAL: TRANSIENTS PLOTS
run ('trp_simu.m')

%% ============================ FUNCTIONS ================================

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

    % Matrix M
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

    % Matrix P
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
    % Forming a symbolic vectors x and x^(2) (x_kr)
    x_G1 = sym('x_g1', [6 1]);
    x_G2 = sym('x_g2', [6 1]);
    x_G3 = sym('x_g3', [6 1]);
    x_G4 = sym('x_g4', [6 1]);
    x = [x_G1; x_G2; x_G3; x_G4];
    x_kr = kron(x,x);
    % Identity matrix
    I = eye(12);
    % Expression for V in symbolic form (see formula (21))
    V = inv(Y - Dd - Sd * P * kron(I,x)) * (Cd * x + Cd_2 * E * kron(x,x));
    % Expressions for the 2nd, 4th and 5th terms of formula (20) in symbolic form
    term_2 = Bd * V;
    term_4 = Bd_2 * M * kron(V,V);
    term_5 = Hd * P * kron(V,x);
    % Decompose term_2, term_4, term_5 into a Taylor series up to the 2nd 
    % order at the point x0 = 0
    TaylorSeries.term2.taylor = taylor(term_2,x,'Order',3);
    TaylorSeries.term4.taylor = taylor(term_4,x,'Order',3);
    TaylorSeries.term5.taylor = taylor(term_5,x,'Order',3);

    % Calculation of matrices from equation (22)

    % Matrices B_D*B_lin & B_D*B_sqr
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
    % Matrices B_D2*M*B_2lin & B_D2*M*B_2sqr
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
    % Matrices H_D*P*H_lin & H_D*P*H_sqr
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
    % Matrices À1, À2
    sys_mat.A1 = Ad + Bd_Blin + Bd2M_B2lin + HdP_Hlin;
    sys_mat.A2 = Ad_2 * E + Bd_Bsqr + Bd2M_B2sqr + HdP_Hsqr;
    % List of states
    sys_mat.x = string(x);
    sys_mat.x_kr = string(x_kr);
    sys_mat.x_full = [sys_mat.x;sys_mat.x_kr];
end

function c_lin = get_lin(c,t,x)
    % get linear coefficients
    c_lin = zeros (1,size(x,1));
    for i_x = 1 : size(x,1)
        id = find(strcmp(x(i_x),t));
        if ~isempty(id)
            c_lin(i_x) = c(id);
        end    
    end
end

function c_sqr = get_sqr (c,t,x_kr)
    % get quadratic coefficients
    exception_x_kr = string;
    c_sqr = zeros (1,size(x_kr,1));
    for i_x = 1 : size(x_kr,1)
        if sum(x_kr(i_x) == exception_x_kr) == 1 % if the element from x_kr has already been
            continue %  then we move on to the next step
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
    % Matrix A of a bilinearized system
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
    % Matrix B (controls) bilinearized system
    B = zeros(size(A,1), 8);
    B(1:n,:) = B1;
    sys_mat.B = B;
    % Matrix N of a bilinearized system
    for i = 1: size(B1,2)
        eval(['B20',num2str(i), ' = kron(B1(:,i), I) + kron(I, B1(:,i));'])
        eval(['sys_mat.N.N',num2str(i), ' = zeros(size(A));'])
        eval(['sys_mat.N.N',num2str(i),'(n+1:n2, 1:n)' ' = B20',num2str(i),';'])
    end
end

function lindep_sys = fun_lindep (sys_mat)
    A1 = sys_mat.A1; % linear part of the bilinear system matrix
    A2 = sys_mat.A2; % quadratic part of the bilinear system matrix 
    B1 = sys_mat.B1;
    N = sys_mat.N;
    % subtract the lines:
    A1([7,8],:) = A1([7,8],:) - A1([1,2],:); 
    A1([13,14],:) = A1([13,14],:) - A1([1,2],:);
    A1([19,20],:) = A1([19,20],:) - A1([1,2],:);
    A1([1,2],:) = A1([1,2],:) - A1([1,2],:);
    A2([7,8],:) = A2([7,8],:) - A2([1,2],:); 
    A2([13,14],:) = A2([13,14],:) - A2([1,2],:);
    A2([19,20],:) = A2([19,20],:) - A2([1,2],:);
    A2([1,2],:) = A2([1,2],:) - A2([1,2],:);
    % reconstructing matrix A21 and matrix A of the bilinearized system
    n = size(A1,1);
    n2 = n + n * n;
    I = eye(n, n);
    A21 = kron(A1, I) + kron(I, A1);
    A_lnd = zeros(n + n^2, n + n^2);
    A_lnd(1:n, 1:n) = A1;
    A_lnd(1:n, n+1:n2) = A2;
    A_lnd(n+1:n2, n+1:n2) = A21;
    % linear and quadratic approximation matrices
    lindep_sys.A1 = A1;
    lindep_sys.A = A_lnd;
    % matrices B & N:
    % subtract the lines:
    B1([7,8],:) = B1([7,8],:) - B1([1,2],:); 
    B1([13,14],:) = B1([13,14],:) - B1([1,2],:);
    B1([19,20],:) = B1([19,20],:) - B1([1,2],:);
    B1([1,2],:) = B1([1,2],:) - B1([1,2],:);
    % Matrix B (controls) bilinearized system
    B = zeros(size(A_lnd,1), 8);
    B(1:n,:) = B1;
    % Matrix N bilinearized system
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
            id(s,1) = p; % list of variable numbers without repetition
        else
            r = r + 1;
            rmv(r,1) = p; % list of duplicate variable numbers to be deleted
            svd(r,1) = tmp(1); % list of duplicate variable numbers to be saved
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