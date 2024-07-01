mu = 0.1;
sigma = 0.2;
s0 = 1;
num = 1000;
long = 1;
rng(123); % 给定一个随机种子，使每次模拟生成的随机数相同
dt = 0.001;

S = zeros(4, num);
S(1,:) = linspace(0, long, num);
S(2:4, 1) = s0; % 初始化S的第二、三、四行的第一个元素

dWt = sqrt(dt) * randn(num, 1);
Wt = cumsum(dWt);

S(2,:) = s0 * exp((mu - 0.5 * sigma^2) * S(1,:) + sigma * Wt.'); % 真实解

for i = 1:num-1
    S(3, i+1) = S(3, i) + mu * S(3, i) * dt + sigma * S(3, i) * dWt(i); % EM方法的解
end

dWt1 = dWt.^2 - dt;
for i = 1:num-1
    S(4, i+1) = S(4, i) + mu * S(4, i) * dt + sigma * S(4, i) * dWt(i) + 0.5 * sigma^2 * S(4, i) * dWt1(i); % Milstein方法的解
end

% 画图对比
% figure(1);
% title('几何布朗运动');
% xlabel('t');
% ylabel('S(t)');
% hold on;
% plot(S(1, 1:1000), S(3, 1:1000), 'r--', 'LineWidth', 1);
% plot(S(1, 1:1000), S(4, 1:1000), 'g--', 'LineWidth', 1);
% plot(S(1, 1:1000), S(2, 1:1000), 'b', 'LineWidth', 1);
% legend('EM', 'MIL', 'True');
% hold off;
% grid on;

% error数组第一行为时间，第二行为em误差，第三行为Mil误差
error = zeros(3, num);
error(1,:) = S(1,:);
error(2,:) = S(3,:) - S(2,:);
error(3,:) = S(4,:) - S(2,:);

% 绘制误差图
figure(2);
xlabel('t');
ylabel('误差');
plot(error(1,:), error(2,:));
hold on
plot(error(1,:), error(3,:));
legend('EM误差', 'MIL误差');


% 计算并输出均方误差
error_EM = mean(error(2,:).^2);
error_MIL = mean(error(3,:).^2);
fprintf('Euler-Maruyama方法的均方误差为：%f\n', error_EM);
fprintf('Milstein方法的均方误差为：%f\n', error_MIL);