
clear all
clear scr
No_cluster_head=5;
pointsNumber=4; 
xx0=[0,1600,-100,1800,1300];
yy0=[0,1900,0,0,-1900] ;
r=50; %radius of disk
R=ones(pointsNumber);
%Simulate binomial point process
theta=2*pi*(rand(pointsNumber,1)); %angular coordinates
rho=r*sqrt(rand(pointsNumber,1)); %radial coordinates
%Convert from polar to Cartesian coordinates
[xx1,yy1]=pol2cart(theta,rho); %x/y coordinates of Poisson points
%Shift centre of disk to (xx0,yy0)
Cluster_head = randi([1 No_cluster_head],1,pointsNumber);
for i=1:pointsNumber
    xx(i)=round(xx1(i)+ xx0(Cluster_head(i)));
    yy(i)=round(yy1(i)+yy0(Cluster_head(i)));
end

%  xx(pointsNumber+1:pointsNumber+No_cluster_head)=xx0;
%  yy(pointsNumber+1:pointsNumber+No_cluster_head)=yy0;
D = [];

for i=1:pointsNumber   
    for j=1:pointsNumber
        Distance_Matrix(i,j)= sqrt((xx(j)-xx(i))^2+(yy(j)-yy(i))^2);
    end
end
% G = digraph(1,1:pointsNumber,D(1));
% for i=2:pointsNumber
%    G = addedge(G,i,1:pointsNumber,D(i));
% end
G = graph(D);
G.Edges.Weight;
number_of_nodes = pointsNumber;
number_of_links = factorial(number_of_nodes - 1);
count = 0;
count_1 = 0;
count_2 =0;
for i=1:pointsNumber
 for j=1:pointsNumber
        D(i,j) = sqrt((xx(j)-xx(i))^2+(yy(j)-yy(i))^2);
 end
end
Recurring_Matrix(1:number_of_nodes) = 0;
%Distance_Threshold = 800;
signal_confirmation=zeros(number_of_nodes);
%Inputing random distances of each link
%Distance_Matrix(1:number_of_nodes,1:number_of_nodes) = randi(1000,9);
%Making diagonal elements of distances zero.
%for i = 1:number_of_nodes
 %   Distance_Matrix(i,i) = 0;
%end
%Making diagonal elements of STP Matrix 1 
% for i = 1:number_of_nodes
%     signal_confirmation(i) = 0;
% end
STP_Matrix =rand(number_of_nodes,number_of_nodes);
%STP_Matrix = fso(Distance_Matrix)
snr=[20];
  for sn=1:length(snr)
count_10 = 0;
 parfor i = 1:number_of_nodes 
     for j = 1:number_of_nodes 
         if Distance_Matrix(i,j)<400
             count_10 = count_10+1;
             STP_Matrix(i,j) = PLC(Distance_Matrix(i,j),snr(sn));
             R(i,j)=0;
         else
             STP_Matrix(i,j) =FSO(Distance_Matrix(i,j),snr(sn));
         end
     
     end
 end
  end
  
 D = zeros(number_of_nodes,number_of_nodes);
 
for i = 1:number_of_nodes 
    STP_Matrix(i,i) = 1;
end

MATRIX = STP_Matrix;
STP_Threshold = 0.70; %Threshold STP
signal_confirmation(1,1)=1;
%Step 1 where the signal is sent to all the nodes connected to node 1, which has STP greater than the threshold value
for i = 2:number_of_nodes
    if STP_Matrix(1,i) > STP_Threshold
       if  signal_confirmation(1,i) == 0
           signal_confirmation(1,i) = 1;
           count = count + 1;
       end 
    end 
end
ASHISH_MC(number_of_nodes,number_of_nodes) = 0;
for i = 2:number_of_nodes
    if signal_confirmation(1,i) == 1
       D(1,i)= sqrt((xx(1)-xx(i))^2+(yy(1)-yy(i))^2);
       D(i,1)= sqrt((xx(i)-xx(1))^2+(yy(i)-yy(1))^2);
       STP_Matrix(i,:) = STP_Matrix(1,i)*STP_Matrix(i,:);
       ASHISH_MC(i,:) = STP_Matrix(i,:);
       for z = 1:number_of_nodes
           if sum(signal_confirmation(:,z)) == 1
               ASHISH_MC(:,z) = 0;
           end
       end
    end
end
SUMMATION = 0;
bhenchod = 0;
con_n = 0;
max_1 = max(max(ASHISH_MC));
if count > 0
    while max_1 > 0.7
        count_1 = count_1 +1;
        for i = 2: number_of_nodes
            for j = 2: number_of_nodes
                if ASHISH_MC(i,j) == max_1
                        if not(i==j)
                            if sum(signal_confirmation(:,j)) == 0 
                            count_2 = count_2 +1;
                            signal_confirmation(i,j) = 1;
                            D(i,j)= sqrt((xx(i)-xx(j))^2+(yy(i)-yy(j))^2);
                            D(j,i)= sqrt((xx(j)-xx(i))^2+(yy(j)-yy(i))^2);
                            STP_Matrix(j,:) = STP_Matrix(i,j)*STP_Matrix(j,:);
                            ASHISH_MC(j,:) = STP_Matrix(j,:);
                            for k = 1:number_of_nodes
                                if sum(signal_confirmation(:,k)) > 0
                                    ASHISH_MC(:,k) = 0;
                                end
                            end
                            max_1 = max(max(ASHISH_MC));
                            else
                            ASHISH_MC(:,j) = 0;
                            end
                        else
                            ASHISH_MC(i,j) = 0;
                        end                    
                end
            end
        end
        max_1 = max(max(ASHISH_MC));
    end
end
%  for g = 1:number_of_nodes
%      if sum(signal_confirmation(:,g)) == 0
%          con_n = con_n + 1; 
%      end
%  end
%  con_ratio(sn)=(pointsNumber-con_n)*100/pointsNumber;
%  con_n=0;
%  end
%  figure(3)
%  plot(snr, con_ratio)
 G = graph(D);
G.Edges.Weight;

p = plot(G)
p.XData = xx;
p.YData = yy;
for i=1:pointsNumber   
    for j=1:pointsNumber
     if R(i,j)==0
         highlight(p,[i,j],'EdgeColor', 'r')
     end
    end
end