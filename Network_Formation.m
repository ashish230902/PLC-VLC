%Network Formation for VLC_PLC Paper
clear all
clear scr

r=100; %radius of disk
xx0=100; yy0=100; %centre of disk
%Simulate binomial point process
pointsNumber=10; 
theta=2*pi*(rand(pointsNumber,1)); %angular coordinates
rho=r*sqrt(rand(pointsNumber,1)); %radial coordinates
%Convert from polar to Cartesian coordinates
[xx,yy]=pol2cart(theta,rho); %x/y coordinates of Poisson points
%Shift centre of disk to (xx0,yy0)
xx=round(xx+xx0);
yy=round(yy+yy0);
%Plotting
D = [];
for i=1:pointsNumber   
 for j=1:pointsNumber
        D(i,j)=sqrt((xx(i)-xx(j))^2+(yy(i)-yy(j))^2);
    end
end
%G = digraph(1,1:pointsNumber,D(1));
%for i=2:pointsNumber
 %   G = addedge(G,i,1:pointsNumber,D(i));
%end
G = graph(D);
G.Edges.Weight;
for i=1:pointsNumber   
    for j=1:pointsNumber
        R(i,j)=randi([0, 1])
    end
end

p = plot(G,'EdgeLabel',G.Edges.Weight)
p.XData = xx;
p.YData = yy;
for i=1:pointsNumber   
    for j=1:pointsNumber
     if R(i,j)==0
         highlight(p, i, j, 'EdgeColor', 'r')
     end
    end
end
