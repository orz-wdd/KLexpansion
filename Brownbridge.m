clc
close all
clear all
N = 20;
% 设置容差和最大迭代次数
tolerance = 1e-5;
maxIter = 1000;
for i = 1:N
    myfun1 = @(x) x + tan(x/2+2*i*pi);
    myfun2 = @(x) 1-x*tan(x/2+2*(i-1)*pi);
    % 选择第K（-1，-2，...个周期）
    a1 = 2*i*pi-1;
    b1 = 2*i*pi-pi+0.01 ;

    a2 = 2*(i-1)*pi;
    b2 = 2*(i-1)*pi+pi-0.01 ;
    % 调用二分法求解
    root1(i) = bisectionMethod(myfun1, a1, b1, tolerance, maxIter);
    root2(i) = bisectionMethod(myfun2, b2, a2, tolerance, maxIter);

    % 输出结果
    %disp(['Root found at:', num2str(root1(i))]);

end
%%
Omega = sort([root1,root2]);

lambda = 2./(1+Omega.^2);
figure(1)

semilogy(lambda(1:10),'o');
xlabel("K-L model index")
ylabel("lambda")

%% 特征函数图像
x_point = 0:0.005:1;
figure(2)
for i = 1:20
    if rem(i, 2) == 1
        ui(i,:) = cos(Omega(i)*(x_point-1/2))/sqrt(1/2+sin(Omega(i))/(2*Omega(i)));
        %disp("222")
    else
        ui(i,:) = sin(Omega(i)*(x_point-1/2))/sqrt(1/2-sin(Omega(i))/(2*Omega(i)));
        %disp("111")
    end
    ylabel("lambda^1/2*ui(x)");
    xlabel("t");
    plot(x_point,sqrt(lambda(i))*ui(i,:));
    hold on

end
N_trunc = 20;
X =zeros(1,length(x_point));
for j = 1:10 %绘制十个样本路径
    rng(j*83)%每次选择不同的样本路径
    for i = 1:N_trunc

        xi = randn;
        X(1,2:201) = X(1,2:201)+xi*sqrt(lambda(i)*ui(i,2:201));
    end
    X(1,1) = 0;
    figure(3)

    plot(x_point,X)
    hold on
end

%%协方差函数的逼近
N_num= [4,8,16];
S = x_point;
T = x_point;
R = zeros(length(S),length(T));
[X,Y] = meshgrid(S,T);
R_ture = exp(-abs(X-Y));
 figure(4)
 subplot(2,2,1)
 mesh(X,Y,R_ture)
 title("真实值")
for j = 1:length(N_num)

    for i = 1:N_num(j)
        R = R+lambda(i)*ui(i,:)'*ui(i,:);
    end
    figure(4)
    subplot(2,2,j+1)

    mesh(X,Y,R)
    titleString = sprintf('截断项为 %d', N_num(j));
    title(titleString)

end

%%
function root = bisectionMethod(f, a, b, tolerance, maxIter)
fa = f(a);
fb = f(b);

if fa * fb > 0
    disp('The chosen interval does not contain a root.');
    return;
end
for i = 1:maxIter
    c = (a + b) / 2;
    fc = f(c);

    if abs(fc) < tolerance
        root = c;
        return;
    end

    if fa * fc < 0
        b = c;
        fb = fc;
    else
        a = c;
        fa = fc;
    end
end
disp('The method did not converge within the maximum number of iterations.');
root = c;
end
