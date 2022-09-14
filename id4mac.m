function [init_data] = id4mac()
% ������ ����� �� Kundur P. "Power Systems Stability and Control". � New 
% York: McGraw-Hill, Inc., 1994. Example 12.6, ���. 813.

% ��������� ����������� � ������� ����� 900MVA/20kV
n_gen = 4; % ���������� �����������
X_d(1:n_gen,1) = 1.80; X_d_ht(1:n_gen,1) = 0.3; X_d_2ht(1:n_gen,1) = 0.25;
X_q(1:n_gen,1) = 1.70; X_q_ht(1:n_gen,1) = 0.55; X_q_2ht(1:n_gen,1) = 0.25;

T_d0_ht(1:n_gen,1) = 8;    T_d0_2ht(1:n_gen,1) = 0.03;
T_q0_ht(1:n_gen,1) = 0.4;  T_q0_2ht(1:n_gen,1) = 0.05;

R_a(1:n_gen,1) = 0.0025;    X_l(1:n_gen,1) = 0.2;

A_Sat(1:n_gen,1) = 0.015; B_Sat(1:n_gen,1) = 9.6; psi_TI(1:n_gen,1) = 0.9;

H = [6.5;6.5;6.175;6.175]; %[G1;G2;G3;G4]

K_D(1:n_gen, 1) = 0;

f = 60; % ������� ���� � ��
omega_0(1:n_gen, 1) = 2*pi*f; % ������� ���� � ��������

init_data.machine = struct('omega_0',omega_0, 'X_d',X_d, 'X_d_ht',X_d_ht,...
    'X_d_2ht',X_d_2ht, 'X_q',X_q, 'X_q_ht',X_q_ht, 'X_q_2ht',X_q_2ht,...
    'T_d0_ht',T_d0_ht,'T_d0_2ht',T_d0_2ht, 'T_q0_ht',T_q0_ht,...
    'T_q0_2ht',T_q0_2ht, 'R_a',R_a, 'X_l',X_l, 'A_Sat',A_Sat,...
    'B_Sat',B_Sat, 'psi_TI',psi_TI, 'H',H, 'K_D',K_D, 'n_gen', n_gen);

% �������� � ���������� ������������ ���:
P = [700;700;719;700]; % �������� �������� � MW, [G1;G2;G3;G4]
Q = [185;235;176;202]; % ���������� �������� � MVAr, [G1;G2;G3;G4]
absE = [1.03;1.01;1.03;1.01]; % ������ ������������ ���������� � pu �������
                              % ����� 100MVA/230kV
angE = [20.2;10.5;-6.8;-17.0];% ���� �������� ������� ���������� � ��������

% �������� ��������:
P_load = [967;1767]; % �������� �������� � MW, [bus7;bus9]
Q_load = [200-100;350-100]; % ���������� �������� � MVAr, [bus7;bus9]
n_load = 2; % ���������� ��� � ���������

init_data.power_flow = struct('P',P, 'Q',Q, 'absE',absE, 'angE',angE,...
    'P_load',P_load, 'Q_load',Q_load, 'n_load',n_load);

% ��������� ����� � ������� �����:
r = 0.0001; % �������� �������� ������������� [p.u. / km]
x_L = 0.001; % �������� ���������� ������������� [p.u / km]
b_C = 0.00175i; % �������� ������������ ������������� ����������
                % ������������� [p.u / km]

r_transf_m = 0.15i; % ���������� ������������� �������������� [p.u.]
                    % � ������� ����� 900MVA/20kV

init_data.line = struct('r',r, 'x_L',x_L, 'b_C',b_C, 'r_transf_m',r_transf_m);
end