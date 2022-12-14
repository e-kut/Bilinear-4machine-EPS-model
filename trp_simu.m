% ????? ???????:
% 1. ???????? ??????????? ???????? ?????????? ????? ?????????? ?????? ?
% ?????????? ????????? ???????????? ??????????

% ????????? ?????????????
dt = 1e-3; % ??? ?????????????
t = 0:dt:50.0; % ????? ?????????????
x0 = zeros(size(sys_mat.init.A1,1),1); % ?????? ?????????? ????????? ??? t = 0
t_size = size(t, 2);


%======== ???? ????????? ?????????? ????? ?????????? ????????? ============
%=(???? ?????????? ???????? ????? ?????????? - ???????????????? ???? ????)=
%==========================================================================
u_1 = zeros(1,length(t));
u_2 = zeros(1,length(t));
x0(8) = 0.5; % ?????? ?????????? ?????????? ?????????
% ???????? ??????? ? ???? ??????? ?????????? ?? ?????? ? ????? ?????????
% ???????????
x0([7,8]) = x0([7,8]) - x0([1,2]);
x0([13,14]) = x0([13,14]) - x0([1,2]);
x0([19,20]) = x0([19,20]) - x0([1,2]);
x0([1,2]) = x0([1,2]) - x0([1,2]);
perturb_index = 1; % ??? ???????????? ????????? ????????
idg = [1;2];
%====== ????? ????? ????????? ?????????? ????? ?????????? ????????? =======
%==========================================================================



% %============== ???? ????????? ?????????? ????? ?????????? ================
% %==(???? ?????????? ???????? ????? ?????? x - ???????????????? ???? ????)==
% %==========================================================================
% % u_1 - ?????????? ???????????? ???????? ?? ?????? ?????????? ? ng
% % u_2 - ?????????? ?????????? ??????? ??????????? ?????????? ? ng
% 
% ng = 2; % ????? ??????????, ????? ???? ???????? ???????? ??????????
% idg = [2*ng-1; 2*ng]; % ?????? "???????????" ?????????? ? ??????? u
% K_u = 1; % ????????? ???????????? ??????????? (??????????)
% %freq_u = 0.0; % ??????? ???????????? ??????????? (??????????) ? ????????
% %u_2 = K_u * sin(freq_u*t); % ?????????? ?????????????? ???????????
% %?????????? ??????? ???????????
% 
% % % ?????????? ??????? ?????????? ??????? ??????????? ?? ???????? K_u ?
% % % ?????????? ??????? ?? t1 ?? t2
% % u_1 = zeros(1,t_size); % ?????????? ???????????? ???????? ?? ?????? ?????????? ? ng
% % u_2 = zeros(1,t_size); % ?????????? ?????????? ??????? ??????????? ?????????? ? ng
% % t1 = 1; % ????? ?????? ?????? ?????????? [?]
% % t2 = 2; % ????? ????????? ?????? ?????????? [?]
% % u_2(t1 / dt : t2 / dt) = K_u;
% % perturb_index = 2; % ??? ???????????? ????????? ????????
% 
% % ?????????? ??????? ???????????? ???????? ?? ???????? K_u ?
% % ?????????? ??????? ?? t1 ?? t2
% u_1 = zeros(1,t_size); % ?????????? ???????????? ???????? ?? ?????? ?????????? ? ng
% t1 = 1; % ????? ?????? ?????? ?????????? [?]
% t2 = 2; % ????? ????????? ?????? ?????????? [?]
% u_1(t1 / dt : t2 / dt) = K_u;
% u_2 = zeros(1,t_size); % ?????????? ?????????? ??????? ??????????? ?????????? ? ng
% perturb_index = 2; % ??? ???????????? ????????? ????????
% %=========== ????? ????? ????????? ?????????? ????? ?????????? ============
% %==========================================================================

%% ???????? ????????:
A_lin = sys_mat.reduce_sys.A1;
B_lin = sys_mat.reduce_sys.B1;
x = x0;
u = zeros (size(B_lin,2),1);
y_lin = zeros(t_size, size(x0,1));

for i = 1:t_size
    y_lin(i, :) = x;
    u(idg(1)) = u_1(i);
    u(idg(2)) = u_2(i);
    dx = A_lin * x + B_lin * u;
    x = x + dx * dt;
end
figure()
plot(t, y_lin)
if perturb_index == 1
    T = '?????????? ??????? ???????? ??????? ??? ?????????? ????? ????????? ???????';
else
    T = '?????????? ??????? ???????? ??????? ??? ?????????? ????? ??????????';
end
title (T)
legend ('\omega_{G1}','\delta_{G1}','\psi_{fd-G1}','\psi_{1d-G1}',...
    '\psi_{1q-G1}','\psi_{2q-G1}','\omega_{G2}','\delta_{G2}',...
    '\psi_{fd-G2}','\psi_{1d-G2}','\psi_{1q-G2}','\psi_{2q-G2}',...
    '\omega_{G3}','\delta_{G3}','\psi_{fd-G3}','\psi_{1d-G3}',...
    '\psi_{1q-G3}','\psi_{2q-G3}','\omega_{G4}','\delta_{G4}',...
    '\psi_{fd-G4}','\psi_{1d-G4}','\psi_{1q-G4}','\psi_{2q-G4}')
xlabel('?????, ?')
%% ?????????? ???????? (???????????? ?????????????):
A_sqr = sys_mat.reduce_sys.A;
B_sqr = sys_mat.reduce_sys.B;
name_field = string(['N' num2str(idg(1))]);
Nu1 = getfield(sys_mat.reduce_sys.N,name_field);
name_field = string(['N' num2str(idg(2))]);
Nu2 = getfield(sys_mat.reduce_sys.N,name_field);
x = [x0; kron(x0,x0)];
u = zeros (size(B_sqr,2),1);
y_bilin2 = zeros(t_size, size(x0,1));
x = x(sys_mat.reduce_sys.ids);

tic
for i = 1:t_size
    y_bilin2(i, :) = x(1:size(x0,1));
    u(idg(1)) = u_1(i);
    u(idg(2)) = u_2(i);
    dx = A_sqr * x + B_sqr * u  + Nu1 * x * u(idg(1)) + Nu2 * x * u(idg(2));
    x = x + dx * dt;
end
toc
figure()
plot(t, y_bilin2)
if perturb_index == 1
    T = '?????????? ??????? ?????????? ??????? ??? ?????????? ????? ????????? ???????';
else
    T = '?????????? ??????? ?????????? ??????? ??? ?????????? ????? ??????????';
end
title (T)
legend ('\omega_{G1}','\delta_{G1}','\psi_{fd-G1}','\psi_{1d-G1}',...
    '\psi_{1q-G1}','\psi_{2q-G1}','\omega_{G2}','\delta_{G2}',...
    '\psi_{fd-G2}','\psi_{1d-G2}','\psi_{1q-G2}','\psi_{2q-G2}',...
    '\omega_{G3}','\delta_{G3}','\psi_{fd-G3}','\psi_{1d-G3}',...
    '\psi_{1q-G3}','\psi_{2q-G3}','\omega_{G4}','\delta_{G4}',...
    '\psi_{fd-G4}','\psi_{1d-G4}','\psi_{1q-G4}','\psi_{2q-G4}')
xlabel('?????, ?')