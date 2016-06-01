% purpose: particle filter for sensor network
%% sensor selection rule : the nearest node by sensor leader be the next
%% sensor leader
% date 2005-11-12

clear all;
clc;
clf;
area=100;%area interested��Ȥ���򳤶�
sensor_number=300;% sensor number�������ڵ���Ŀ
rand('state',0);
randn('state',0);
dis_sensor=area*rand(sensor_number,2);%distrubtion of sensor nodes
figure(1)
plot(dis_sensor(:,1),dis_sensor(:,2),'bo');

x0=[15,25];%inti of target position��ʼĿ��λ��

%%%%%%%%% start of the algorithm------�㷨

%inti��ʼ��

sigmaT =0.2^2;              % Variance of the Gaussian transition prior.state eqution noise  ת��״̬�������ֵ��˹������
sigmaM =0.1^2;                % Variance of the Gaussian likelihood.measurement noise��˹��Ȼ�۲�����
N = 100;                    % Number of particles.������Ŀ
T =50;                      % Number of time steps.say.the number of iteration����ʱ��

%plot the figure of target ��Ŀ��
v0=zeros(T+1,2);
v0(1,:)=[1.2,0.3];%inti the velocity��ʼ��Ŀ���ٶ�
x_real=zeros(T+1,2);
x_real(1,:)=x0;%Ŀ���һ��������λ��
for t=2:T+1                  % real trajectory of targetʵ��Ŀ��켣
    v0(t,1)=v0(t-1,1)+sqrt(sigmaT)*randn(1);
    v0(t,2)=v0(t-1,2)+sqrt(sigmaT)*randn(1);
    x_real(t,:)=x_real(t-1,:)+v0(t,:);
end

%distributed particle filter�ֲ�ʽ�����˲�

x = zeros(T+1,N,2);%��ͬʱ��N�����ӵ�λ��
x_mmse=zeros(T,2); %��ͬʱ��N�����ӵ�����
velocity=zeros(T+1,N,2);%��ͬʱ��N�����ӵ��ٶ�
vrandn=sqrt(sigmaT)*randn(1,N);%��˹����ֲ�
velocity(1,:,1)=v0(1,1)*ones(1,N)+vrandn;
velocity(1,:,2)=v0(1,2)*ones(1,N)+vrandn;
z = zeros(T,1);%�۲�
u = zeros(T,N,2);%processing noise��������
v = zeros(T,N);%measure noise�۲�����
xrandn=randn(1,N);
x(1,:,1)=x0(1,1)*ones(1,N)+xrandn;%prior dist
x(1,:,2)=x0(1,2)*ones(1,N)+xrandn;
w=zeros(T+1,N);
w(1,:)=rand(1,N);%prior dist
w(1,:)=w(1,:)./sum(w(1,:));%Ȩ�ع�һ��
%inti of the leader the nearest node of target ��ʼ������Ŀ��������ű�ڵ�
L=length(dis_sensor(:,1));
sensor_leader=zeros(T,1);%ͷ������
for k=1:L
    dist1(k)=sqrt((dis_sensor(k,1)-x0(1,1))^2+(dis_sensor(k,2)-x0(1,2))^2);
end
[ss,sensor_leader(1)]=min(dist1);
sensor_range=10;%sensor range ��֪�뾶
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:T
  u(t,:,1)=sqrt(sigmaT)*randn(1,N);                      % Process noise.��������
  u(t,:,2)=sqrt(sigmaT)*randn(1,N);
  velocity(t,:,:)=velocity(t-1,:,:)+u(t,:,:);%%sample from the important function
  x(t,:,:)=x(t-1,:,:)+velocity(t,:,:);
    
  %select sensor leader based on NN rule ����NN����ѡ���ű�ڵ�
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
  plot(dis_sensor(sensor_leader(t),1),dis_sensor(sensor_leader(t),2),'k+');%plot the inti sensor leader ����ʼ�ű�ڵ�
  hold off;
    
  %%% important function update   ��Ҫ�Ժ�������
  z(t,1)=atan((x_real(t,2)-dis_sensor(sensor_leader(t),2))/((x_real(t,1)-dis_sensor(sensor_leader(t),1))))+sqrt(sigmaM)*randn;%leader sensor measurement
  for jj=1:N
      w(t,jj)=w(t-1,jj)*exp(-1/(2*sigmaM)*(z(t,1)-atan((x(t,jj,2)-dis_sensor(sensor_leader(t),2))/((x(t,jj,1)-dis_sensor(sensor_leader(t),1)))))^2);
  end
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.Ȩ�ع�һ��

  %resampling �ز���
  % SELECTION STEP: ѡ�񲽳�
  
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
