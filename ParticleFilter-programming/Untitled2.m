% [V,C]=VoronoiLimit(x,y)
% Provides the Voronoi decomposition of a set of data, but with all
% vertices limited to the boundary created by the data itself. 
% V contains all vertices and C contains all vertices for each individual
% point. That is: V(C{ij},:) will give you the vertices of the ij'th data point. 
%
% Requires the Polybool function of the mapping toolbox to run. 
%
% Run with no input to see example
%
% Made by: Jakob Sievers, Jakob.Sievers@gmail.com
clc;
clear all;
tic
%% 输入语句

   val=5;
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
    C{ij}
end
%%
figure(1)
hold on
plot(x,y,'.k')
voronoi(x,y);
xzuobi=[0,100,100,0,0];
yzuobi=[0,0,100,100,0];
plot(xzuobi,yzuobi,'-b');
axis([bnd(1) bnd(2) bnd(3) bnd(4)+0.1])
xlabel('x/m');
ylabel('y/m');
title('Limited Voronoi Decomposition')
%%
figure(2)
hold on
plot(x,y,'.k')
voronoi(x,y);
for id=1:length(C)
    plot(V(C{id},1),V(C{id},2),'-b')

  numv = size(V,1);
  vlabels = arrayfun(@(n) {sprintf('V%d', n)}, (2:numv)');
  hold on
  Hpl = text(V(2:end,1), V(2:end,2)+.2, vlabels, ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
  grid on;
end

axis([bnd(1) bnd(2) bnd(3) bnd(4)+0.1])
xlabel('x/m');
ylabel('y/m');
title('Limited Voronoi Decomposition (Blue)')
% dx=(bnd(2)-bnd(1))/10;
% dy=(bnd(4)-bnd(3))/10;
% axis([bnd(1)-dx bnd(2)+dx bnd(3)-dy bnd(4)+dy])
% title('Original Voronoi Decomposition (Blue), vs. New limited Voronoi Decomposition (Red)')

nump = size(x,1);
plabels = arrayfun(@(n) {sprintf('X%d', n)}, (1:nump)');
hold on
Hp2 = text(x, y+0.01, plabels, 'color', 'r', ...
      'FontWeight', 'bold', 'HorizontalAlignment',...
      'center', 'BackgroundColor', 'none');
% % Rs=30;
% % w=0:pi/50:2*pi;
% %   x11=x(15)+Rs*cos(w);
% %     y11=y(15)+Rs*sin(w);%表示圆每一个角度对应一个点，然后把这些点连接起来构成圆
% %     hold on
%     plot(x11,y11,'r');   %画出各个传感器节点的感知范围
