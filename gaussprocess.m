clc
close all
clear all
N = 20;
% 设置容差和最大迭代次数
tolerance = 1e-5;
maxIter = 1000;
x_point = 0:0.005:1;
for i = 1:N
    myfun1 = @(x) sqrt(2)*sin(i*pi*x);
    lambda(i) = 1/((pi*i)^2);
    ui(i,:) = myfun1(x_point);

end

figure(1)

semilogy(lambda(1:10),'o');
xlabel("K-L model index")
ylabel("lambda")

%% 特征函数图像

figure(2)
for i = 1:20
  
       
      
    
    ylabel("sqrt(\lambda)*ui(x)");
    xlabel("t");
    plot(x_point,sqrt(lambda(i))*ui(i,:));
    hold on

end
%% 
N_trunc = 20;
X =zeros(1,length(x_point));
for j = 1:10 %绘制十个样本路径
    rng(j)%每次选择不同的样本路径
    for i = 1:N_trunc

        xi = randn;
        X(1,2:200) = X(1,2:200)+xi*sqrt(lambda(i))*ui(i,2:200);
    end
    X(1,1) = 0;
    X(1,end) = 0;
    figure(3)

    plot(x_point,X)
    hold on
end

%% 协方差函数的逼近
N_num= [4,8,16];
S = x_point;
T = x_point;
R = zeros(length(S),length(T));
[X,Y] = meshgrid(S,T);
R_ture = min(X,Y)-X.*Y;
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

