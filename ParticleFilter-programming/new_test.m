
clc;
clear all;
tic
%% 输入语句
   val=300;
   x=100*rand(val,1);
   y=100*rand(val,1);
   xy=unique([x,y],'rows');
   x=xy(:,1);
   y=xy(:,2);


if ~any(size(x)==1) || ~any(size(y)==1) || numel(x)==1 || numel(y)==1
    disp('Input vectors should be single rows or columns')
    return
end

x=x(:);
y=y(:);
bnd=[0 100 0 100]; 
% bnd=[min(x) max(x) min(y) max(y)]; %data bounds
crs=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]); %data boundary corners
%% 划分
dt=DelaunayTri(x(:),y(:));
[V,C]=voronoiDiagram(dt); %This structure gives vertices for each individual point but is missing all "infinite" vertices
[vx,vy]=voronoi(x,y); %This structure includes the "infinite" vertices but provides everything as a completele list of vertices rather than individually for each point. 
%Hence we need to add the missing vertices from vx and vy to the V and C structure. 
vxyl=[vx(:) vy(:)];
xix=ones(size(vx));

%Eliminate spurious double entries in V(C{ij})
av=cell(length(C),1);
for ih=1:length(C)
    av{ih}=[V(C{ih},1),V(C{ih},2)];  
    if size(unique(av{ih},'rows'),1)<size(av{ih},1);
        k=0;
        ctr=0;
        while k==0
            ctr=ctr+1;
            if ctr==1
                tt=unique(V(C{ih}(ctr),:)==V(C{ih}(end),:));
                if length(tt)<2 && tt==1
                    C{ih}(end)=[];
                    k=1;
                else
                    tt=unique(V(C{ih}(ctr),:)==V(C{ih}(ctr+1),:));
                    if length(tt)<2 && tt==1
                        C{ih}(ctr+1)=[];
                        k=1;
                    end
                end
            else
                tt=unique(V(C{ih}(ctr),:)==V(C{ih}(ctr+1),:));
                if length(tt)<2 && tt==1
                    C{ih}(ctr+1)=[];
                    k=1;
                end
            end
        end
    end
end
lV0=length(V);

%Find missing points that should be added to existing V/C structure
for ii=1:length(vxyl)
    fix=find(V(:,1)==vxyl(ii,1));
    if ~isempty(fix)
        if any(V(fix,2)==vxyl(ii,2))
            xix(ii)=0;
        end
    end
end

mix=find(xix==1)./2; %index of missing values
lmix=length(mix);
mvx=vx(2,mix); %missing vx
mvy=vy(2,mix); %missing vy
mv=[mvx',mvy'];
cpx=vx(1,mix); %connector point x (connects between outer missing points and inner existing points in V/C)
cpy=vy(1,mix); %connector point y (connects between outer missing points and inner existing points in V/C)

ctr=0;
mv2=[];
cpVixt=cell(lmix,1); %connector points, index in V structure
for ii=1:lmix
    if any(V(:,1)==cpx(ii) & V(:,2)==cpy(ii))
        cpVixt{ii}=find(V(:,1)==cpx(ii) & V(:,2)==cpy(ii));
        lval=length(cpVixt{ii});
        if lval==1
            ctr=ctr+1;
            mv2(ctr,:)=mv(ii,:);
        elseif lval>1
            ctr=ctr+1;
            mv2(ctr:ctr+lval-1,:)=[ones(lval,1).*mv(ii,1) ones(lval,1).*mv(ii,2)];
            ctr=ctr+lval-1;
        end
    end
end
cpVixt=cell2mat(cpVixt);


V=[V;mv2]; %add points to V structure
allVixinp=inpolygon(V(:,1),V(:,2),crs(:,1),crs(:,2)); %determine which points in V that are within the data boundaries. 

%Addition-routine: addition of missing points (mvx,mvy) to individual vertice-polygons (C)
for ij=1:length(C)
    if any(C{ij}==1)
        ixa=find(cpVixt==C{ij}(2));
        ixb=find(cpVixt==C{ij}(end));
        if length(C{ij})<3 %corner point detected
            C{ij}(1)=lV0+ixa(1);
            C{ij}=[C{ij},lV0+ixa(2)];
        else
            if length(ixa)==1 && length(ixb)==1
                C{ij}(1)=lV0+ixa;
                C{ij}=[C{ij},lV0+ixb];
            elseif length(ixa)==2 && length(ixb)==1
                C{ij}=[C{ij},lV0+ixb];
                [~,minix]=min(sqrt((V(C{ij}(end),1)-V(cpVixt(ixa),1)).^2+(V(C{ij}(end),2)-V(cpVixt(ixa),2)).^2));
                C{ij}(1)=lV0+ixa(minix);
            elseif length(ixa)==1 && length(ixb)==2
                C{ij}(1)=lV0+ixa;
                [~,minix]=min(sqrt((V(C{ij}(1),1)-V(cpVixt(ixb),1)).^2+(V(C{ij}(1),2)-V(cpVixt(ixb),2)).^2));
                C{ij}=[C{ij},lV0+ixb(minix)];
            elseif length(ixa)==2 && length(ixb)==2
                [~,minix1]=min(sqrt((x(ij)-V(lV0+ixa,1)).^2+(y(ij)-V(lV0+ixa,2)).^2));
                [~,minix2]=min(sqrt((x(ij)-V(lV0+ixb,1)).^2+(y(ij)-V(lV0+ixb,2)).^2));
                C{ij}(1)=lV0+ixa(minix1);
                C{ij}=[C{ij},lV0+ixb(minix2)];
            end
        end
    end
end
%Polybool for restriction of polygons to domain.
C1=C; %Do this analysis based on old vertice descriptions to avoid problems
for ij=1:length(C)
    if sum(allVixinp(C{ij}))~=length(C{ij})
        [xb, yb] = polybool('intersection',crs(:,1),crs(:,2),V(C1{ij},1),V(C1{ij},2));
        ix=nan(1,length(xb));
        for il=1:length(xb)
            if any(V(:,1)==xb(il)) && any(V(:,2)==yb(il))
                ix1=find(V(:,1)==xb(il));
                ix2=find(V(:,2)==yb(il));
                for ib=1:length(ix1)
                    if any(ix1(ib)==ix2)
                        ix(il)=ix1(ib);
                    end
                end
                if isnan(ix(il))==1
                    lv=length(V);
                    V(lv+1,1)=xb(il);
                    V(lv+1,2)=yb(il);
                    allVixinp(lv+1)=1;
                    ix(il)=lv+1;
                end
            else
                lv=length(V);
                V(lv+1,1)=xb(il);
                V(lv+1,2)=yb(il);
                allVixinp(lv+1)=1;
                ix(il)=lv+1;
            end
        end
        C{ij}=ix;
    end
    C{ij};
end
%%
% figure(1)
% hold on
% voronoi(x,y);
% plot(x,y,'co')
% xzuobi=[0,100,100,0,0];
% yzuobi=[0,0,100,100,0];
% plot(xzuobi,yzuobi,'-b');
% axis([bnd(1) bnd(2) bnd(3) bnd(4)+0.1])
% xlabel('x/m');
% ylabel('y/m');
% title('Limited Voronoi Decomposition')


%%
figure(2)
hold on
plot(x,y,'c.');
%%%%
Rs=20;
w=0:pi/50:2*pi;
fugai=zeros(val,1);
flag=zeros(val,1);
num=0;
xx0=[];
for i=1:val
    for j=1:val
        if i~=j && flag(i)==0
            dis=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
            if dis<=Rs
                dis1=sqrt((x(i)-V(C{j},1)).^2+(y(i)-V(C{j},2)).^2);
                n=max(dis1);
                if n<=Rs
                    fugai(j)=1;
                    flag(j)=1;
                end                
            end
        end
     end
     if fugai(i)==0 && flag(i)==0
         x11=x(i)+Rs*cos(w);
         y11=y(i)+Rs*sin(w);
         hold on
%          plot(x11,y11,'r:'); 
         num=num+1;
         xx0=[xx0;x(i),y(i)];
     end
end
num
xx0;
 voronoi(xx0(:,1),xx0(:,2),'go')
  plot(xx0(:,1),xx0(:,2),'bo')
    
xzuobi=[0,100,100,0,0];
yzuobi=[0,0,100,100,0];
plot(xzuobi,yzuobi,'-b');
axis equal;
axis([bnd(1) bnd(2) bnd(3) bnd(4)+0.1])
xlabel('x/m');
ylabel('y/m');
title('Voronoi Coverage')




% purpose: particle filter for sensor network
%% sensor selection rule : the nearest node by sensor leader be the next
%% sensor leader
% date 2005-11-12

% clear all;
% clc;
% clf;
% area=100;%area interested兴趣区域长度
% sensor_number=300;% sensor number传感器节点数目
rand('state',0);
randn('state',0);
% dis_sensor=area*rand(sensor_number,2);%distrubtion of sensor nodes

dis_sensor(:,1)=x(:);
dis_sensor(:,2)=y(:);

dis_sensor1=xx0(:,:);



figure(3)

plot(dis_sensor(:,1),dis_sensor(:,2),'bo');
axis equal;
hold on

x0=[15,25];%inti of target position初始目标位置

%%%%%%%%% start of the algorithm------算法

%inti初始化

sigmaT =0.05^2;              % Variance of the Gaussian transition prior.state eqution noise  转移状态方程零均值高斯白噪音
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
u = zeros(T,N,2);%processing noise系统噪音
v = zeros(T,N);%measure noise观测噪音
xrandn=randn(1,N);
x(1,:,1)=x0(1,1)*ones(1,N)+xrandn;%prior dist
x(1,:,2)=x0(1,2)*ones(1,N)+xrandn;
w=zeros(T+1,N);
w(1,:)=rand(1,N);%prior dist
w(1,:)=w(1,:)./sum(w(1,:));%权重归一化
%inti of the leader the nearest node of target 初始化距离目标最近的信标节点
L=length(dis_sensor(:,1));
sensor_leader=zeros(T,1);
for k=1:L
    dist1(k)=sqrt((dis_sensor(k,1)-x0(1,1))^2+(dis_sensor(k,2)-x0(1,2))^2);
end
[ss,sensor_leader(1)]=min(dist1);
sensor_range=10;%sensor range 感知半径
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:T
  
  u(t,:,1)=sqrt(sigmaT)*randn(1,N);                      % Process noise.
  u(t,:,2)=sqrt(sigmaT)*randn(1,N);
  velocity(t,:,:)=velocity(t-1,:,:)+u(t,:,:);%%sample from the important function
  x(t,:,:)=x(t-1,:,:)+velocity(t,:,:)+u(t,:,:);
    
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
  if ss>sensor_range
      disp('The range of sensor is so small, please enlarge the sensor range!  The min sensor range is');
      disp(num2str(ss));
  end
      
  hold on;
%   t
%   pause(0.02);
  
%   plot(dis_sensor(sensor_leader(t),1),dis_sensor(sensor_leader(t),2),'k+');%plot the inti sensor leader 画初始信标节点
%   hold off;
  
  %%% important function update
  z(t,1)=atan((x_real(t,2)-dis_sensor(sensor_leader(t),2))/((x_real(t,1)-dis_sensor(sensor_leader(t),1))))+sqrt(sigmaM)*randn;%leader sensor measurement
  for jj=1:N
      w(t,jj)=w(t-1,jj)*exp(-1/(2*sigmaM)*(z(t,1)-atan((x(t,jj,2)-dis_sensor(sensor_leader(t),2))/((x(t,jj,1)-dis_sensor(sensor_leader(t),1)))))^2);
  end
  w(t,:) = w(t,:)./sum(w(t,:));                % Normalise the weights.

  %resampling 重采样
  % SELECTION STEP:
  
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
figure(3);
hold on;
axis equal;
axis([bnd(1) bnd(2) bnd(3) bnd(4)+0.1])
xlabel('x/m');
ylabel('y/m');
plot(x_real(:,1),x_real(:,2),'b-',x_mmse(2:T,1),x_mmse(2:T,2),'r:');%plot the true trajectory  %plot the true trajectory
% legend('sensor','leader node','Actual trajectory','Estimated trajectory');
% xlabel('x','fontsize',15);
% ylabel('y','fontsize',15);
title('performance of tracking target','fontsize',15);
   
% xzuobi=[0,100,100,0,0];
% yzuobi=[0,0,100,100,0];
% plot(xzuobi,yzuobi,'-b');


figure(4);
hold on;
plot(1:50,error_(:));
xlabel('time')
ylabel('error')

%%
toc


