x = 0; % 初始状态?
R = input('请输入过程噪声方差R的值:'); % 测量噪声协方差
Q = input('请输入观测噪声方差Q的值: '); % 过程状态协方差
tf = 100; % 模拟长度
N = 100; % 粒子滤波器中的粒子个数
xhat = x;
P = 2;
xhatPart = x;

%粒子滤波器初期
for i = 1 : N
    xpart(i) = x + sqrt(P) * randn;
end
xArr = [x];
yArr = [x^2 / 20 + sqrt(R) * randn];
xhatArr = [x];
PArr = [P];
xhatPartArr = [xhatPart];
close all;

for k = 1 : tf
    % 模拟系统
    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    y = x^2 / 20 + sqrt(R) * randn;
    % 粒子滤波器
    for i = 1 : N
        xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
        ypart = xpartminus(i)^2 / 20;
        vhat = y - ypart;% 观测和预测的差?
        q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
    end
    % 平均每一个估计的可能性
    qsum = sum(q);
    for i = 1 : N
        q(i) = q(i) / qsum;%归一化权值
    end
    % 重采样
    for i = 1 : N
        u = rand;
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + q(j);
            if qtempsum >= u
                xpart(i) = xpartminus(j);
                break;
            end
        end
    end
   
    xhatPart = mean(xpart);
   
   % 在图片中表示出数据
    xArr = [xArr x];
    yArr = [yArr y];
    xhatArr = [xhatArr xhat];
    PArr = [PArr P];
    xhatPartArr = [xhatPartArr xhatPart];
end

t = 0 : tf;

figure;
plot(t, xArr, 'b.', t, xhatPartArr, 'k-');
% set(gca,'FontSize',10); set(gcf,'Color','yellow');
xlabel('time step'); ylabel('state');
legend('真实值','PF估计值');

xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
disp(['PF估计误差均方值 = ', num2str(xhatPartRMS)]);
