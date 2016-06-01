% purpose: particle filter for sensor network
%% sensor selection rule : the nearest node by sensor leader be the next
%% sensor leader
% date 2005-11-12

clear all;
clc;
clf;
area=100;%area interested兴趣区域长度
sensor_number=300;% sensor number传感器节点数目
rand('state',0);
randn('state',0);
dis_sensor=area*rand(sensor_number,2);%distrubtion of sensor nodes
figure(1)
plot(dis_sensor(:,1),dis_sensor(:,2),'bo');

x0=[15,25];%inti of target position初始目标位置

%%%%%%%%% start of the algorithm------算法

%inti初始化

sigmaT =0.2^2;              % Variance of the Gaussian transition prior.state eqution noise  转移状态方程零均值高斯白噪音
sigmaM =0.1^2;                % Variance of the Gaussian likelihood.measurement noise高斯似然观测噪音
N = 100;                    % Number of particles.粒子数目
T =50;                      % Number of time steps.say.the number of iteration迭代时间

%plot the figure of target 画目标
v0=zeros(T+1,2);
v0(1,:)=[1.2,0.3];%inti the velocity初始化目标速度
x_real=zeros(T+1,2);
x_real(1,:)=x0;%目标第一个采样点位置
for t=2:T+1                  % real trajectory of target实际目标轨迹
    v0(t,1)=v0(t-1,1)+sqrt(sigmaT)*randn(1);
    v0(t,2)=v0(t-1,2)+sqrt(sigmaT)*randn(1);
    x_real(t,:)=x_real(t-1,:)+v0(t,:);
end

%distributed particle filter分布式粒子滤波

x = zeros(T+1,N,2);%不同时刻N个粒子的位置
x_mmse=zeros(T,2); %不同时刻N个粒子的噪音
velocity=zeros(T+1,N,2);%不同时刻N个粒子的速度
vrandn=sqrt(sigmaT)*randn(1,N);%高斯随机分布
velocity(1,:,1)=v0(1,1)*ones(1,N)+vrandn;
velocity(1,:,2)=v0(1,2)*ones(1,N)+vrandn;
z = zeros(T,1);%观测
u = zeros(T,N,2);%processing noise生成噪音
v = zeros(T,N);%measure noise观测噪音
xrandn=randn(1,N);
x(1,:,1)=x0(1,1)*ones(1,N)+xrandn;%prior dist
x(1,:,2)=x0(1,2)*ones(1,N)+xrandn;
w=zeros(T+1,N);
w(1,:)=rand(1,N);%prior dist
w(1,:)=w(1,:)./sum(w(1,:));%权重归一化
%inti of the leader the nearest node of target 初始化距离目标最近的信标节点
L=length(dis_sensor(:,1));
sensor_leader=zeros(T,1);%头结点矩阵
for k=1:L
    dist1(k)=sqrt((dis_sensor(k,1)-x0(1,1))^2+(dis_sensor(k,2)-x0(1,2))^2);
end
[ss,sensor_leader(1)]=min(dist1);
sensor_range=10;%sensor range 感知半径
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:T
  u(t,:,1)=sqrt(sigmaT)*randn(1,N);                      % Process noise.生成噪音
  u(t,:,2)=sqrt(sigmaT)*randn(1,N);
  velocity(t,:,:)=velocity(t-1,:,:)+u(t,:,:);%%sample from the important function
  x(t,:,:)=x(t-1,:,:)+velocity(t,:,:);
    
  %select sensor leader based on NN rule 根据NN规则选择信标节点
  for k=1:L
      for kk=1:t-1%be sure that the sensor can't be use as leader again
          if ~(k==sensor_leader(kk))
            dist2(k)=sqrt((dis_sensor(k,1)-dis_sensor(sensor_leader(t-1),1))^2+(dis_sensor(k,2)-dis_sensor(sensor_leader(t-1),2))^2);
          end
          dist2(sensor_leader(kk))=10000;    
      end
  end

  [ss,sensor_leader(t)]=min(dist2);
%   if ss>sensor_range
%       disp('The range of sensor is so small, please enlarge the sensor range!  The min sensor range is');
%       disp(num2str(ss));
%   end
      
  hold on;
  t
  pause(0.01);
  plot(dis_sensor(sensor_leader(t),1),dis_sensor(sensor_leader(t),2),'k+');%plot the inti sensor leader 画初始信标节点
  hold off;
    
  %%% important function update   重要性函数更新
  z(t,1)=atan((x_real(t,2)-dis_sensor(sensor_leader(t),2))/((x_real(t,1)-dis_sensor(sensor_leader(t),1))))+sqrt(sigmaM)*randn;%leader sensor measurement
  for jj=1:N
      w(t,jj)=w(t-1,jj)*exp(-1/(2*sigmaM)*(z(t,1)-atan((x(t,jj,2)-dis_sensor(sensor_leader(t),2))/((x(t,jj,1)-dis_sensor(sensor_leader(t),1)))))^2);
  end
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.权重归一化

  %resampling 重采样
  % SELECTION STEP: 选择步长
  
  Nth=N/3;
  Neff=1/(sum(w(t,:).^2)+eps);
  if Neff<Nth               %%resampling
         c=zeros(N,1);
         c(1,1)=w(t,1);
         for i=2:N
             c(i,1)=c(i-1,1)+w(t,i);     
         end
         i=1;
         uu=zeros(N,1);
         uu(1,1)=1/N*rand;
         outIndex=zeros(1,N);
         for j=1:N
             uu(j,1)=uu(1,1)+1/N*(j-1);
             while uu(j,1)>c(i)
                 i=i+1;
             end
             x(t,j,:)=x(t,i,:);
             w(t,j)=1/N;
         end
  end
end

%%plot the result 
                     
for t=2:T
    x_mmse(t,1)=sum(w(t,:).*x(t,:,1));
    x_mmse(t,2)=sum(w(t,:).*x(t,:,2));
    for j=1:N
        error1(t,j)=w(t,j)*((x(t,j,1)-x_real(t,1))^2+(x(t,j,2)-x_real(t,2))^2);
    end
    error_(t)=sqrt(sum(error1(t,:)));
end

error_mean=1/T*sum(error_(:));
figure(1);
hold on;
plot(x_real(:,1),x_real(:,2),'b-',x_mmse(2:T,1),x_mmse(2:T,2),'r:');%plot the true trajectory  %plot the true trajectory
% legend('sensor','leader node','Actual trajectory','Estimated trajectory');
xlabel('x','fontsize',15);
ylabel('y','fontsize',15);
title('performance of tracking target','fontsize',15);
hold off;      

figure(2);
hold on;
plot(1:50,error_(:));
xlabel('time')
ylabel('error')
