x = 0; % ��ʼ״̬?
R = input('�����������������R��ֵ:'); % ��������Э����
Q = input('������۲���������Q��ֵ: '); % ����״̬Э����
tf = 100; % ģ�ⳤ��
N = 100; % �����˲����е����Ӹ���
xhat = x;
P = 2;
xhatPart = x;

%�����˲�������
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
    % ģ��ϵͳ
    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    y = x^2 / 20 + sqrt(R) * randn;
    % �����˲���
    for i = 1 : N
        xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
        ypart = xpartminus(i)^2 / 20;
        vhat = y - ypart;% �۲��Ԥ��Ĳ�?
        q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
    end
    % ƽ��ÿһ�����ƵĿ�����
    qsum = sum(q);
    for i = 1 : N
        q(i) = q(i) / qsum;%��һ��Ȩֵ
    end
    % �ز���
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
   
   % ��ͼƬ�б�ʾ������
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
legend('��ʵֵ','PF����ֵ');

xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
disp(['PF����������ֵ = ', num2str(xhatPartRMS)]);
