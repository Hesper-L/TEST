%本函数利用二分法求解X = ls(Xk,Pk)问题

%目标函数：f

%符号参数：var

%终止限：eps

function x = dichotomy(f,var,eps)

g = diff(f,var);

[a, b] = search(f,var);

x = (a + b)/2; %防止eps过大导致x无值

while b - a > eps

x = (a+b)/2;

gx = subs(g, var, x);

if gx > 0

b = x;

elseif gx < 0

a = x;

else break;

end

end

%加步搜索法-确定搜索区间

function [a, b] = search(g,var)

gt = matlabFunction(g);

X = 0; tmp = X;

h = 1; k = 0;

while 1

Xk = X + h; k = k+1;

Y = subs(gt,var,X);

Yk = subs(gt,var,Xk);

if Y > Yk %加大步长搜索

h = 2 * h;

tmp = X;

X = Xk;

elseif Y == Yk %缩小步长搜索

h = h/2;

elseif k == 1

h = -h; %反向搜索

else break;

end

end

a = min(tmp, Xk);

b = max(tmp, Xk);

end

end

%DFP算法的实现

%有了一维搜索函数，那么实现DFP算法也就能依照算法流程图来设计了：

%DFP算法主程序

%目标函数：f

%初始点：X0

%参数：var

%终止限：eps

function DFP(f, X0, var, eps)

%初始化符号函数，梯度，维数等

syms var t;

g = jacobian(f)'; %Jacobian转置->Grad

fx = matlabFunction(f); %符号函数->函数句柄(R2009以上支持)

gx = matlabFunction(g);

n = length(var); %维数

X = X0; Xk = X0;

while 1

fx0 = fx(X(1),X(2)); gx0 = gx(X(1),X(2));

Hk = eye(2); Pk = -gx0; %初始方向

k = 0; %迭代次数

while 1

Y = Xk + t*Pk;

y = fx(Y(1),Y(2));

tk = dichotomy(y, t, eps); %一维搜索

Xk = Xk + tk*Pk;

fx1 = fx(Xk(1),Xk(2));

gx1 = gx(Xk(1),Xk(2));

if norm(gx1) < eps || k == n

X = Xk; fx0 = fx1;

break;

end

Sk = Xk - X; Yk = gx1 - gx0;

Hk = Hk + Sk*Sk'/(Sk'*Yk) - Hk*(Yk)*Yk'*Hk/(Yk'*Hk*Yk);

Pk = -Hk*gx1; %校正方向

k = k+1;

end

if norm(gx1) < eps

disp('X(k+1) = '); disp(Xk);

disp('F(K+1) = '); disp(fx0);

break;

end