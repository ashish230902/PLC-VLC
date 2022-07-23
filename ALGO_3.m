%Network Formation for VLC_PLC Paper
clear all
clear scr
r=2000; %radius of disk
xx0=100; yy0=100; %centre of disk
%Simulate binomial point process
pointsNumber=6; 
theta=2*pi*(rand(pointsNumber,1)); %angular coordinates
rho=r*sqrt(rand(pointsNumber,1)); %radial coordinates
%Convert from polar to Cartesian coordinates
[xx,yy]=pol2cart(theta,rho); %x/y coordinates of Poisson points
%Shift centre of disk to (xx0,yy0)
xx=round(xx+xx0);
yy=round(yy+yy0);
xx=[0,-25,0,-1000,0,-1000];
yy=[0,1500,1500,0,1700,1500];

%Plotting
D = [];
R=ones(pointsNumber);

for i=1:pointsNumber   
    for j=1:pointsNumber
        D(i,j)= 0;
    end
end
%G = digraph(1,1:pointsNumber,D(1));
%for i=2:pointsNumber
 %   G = addedge(G,i,1:pointsNumber,D(i));
%end
G = graph(D);
G.Edges.Weight;
number_of_nodes = pointsNumber;
number_of_links = factorial(number_of_nodes - 1);
count = 0;
count_1 = 0;
count_2 =0;
for i = 1:number_of_nodes
    for j = 1:number_of_nodes
        Distance_Matrix(i,j) = sqrt((xx(j)-xx(i))^2+(yy(j)-yy(i))^2);
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

% 
 for i = 1:number_of_nodes 
     for j = 1:number_of_nodes 
         if Distance_Matrix(i,j)<400
             STP_Matrix(i,j) = PLC(Distance_Matrix(i,j));
             R(i,j)=0;
         else
             STP_Matrix(i,j) =FSO(Distance_Matrix(i,j));
         end
             
     
     end
 end
 
 
for i = 1:number_of_nodes 
    STP_Matrix(i,i) = 0;
end
MATRIX = STP_Matrix;
STP_Threshold = 0.7; %Threshold STP
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
       ASHISH_MC_2 = ASHISH_MC;
    end
end
blah(number_of_nodes) = 0;
why = 9;
%STEP 2 

max_1 = max(max(ASHISH_MC));
while max_1 > STP_Threshold
    count_1 = count_1 +1;
    for i = 2: number_of_nodes
        for j = 2: number_of_nodes
            if ASHISH_MC(i,j) == max_1
                if not(i==j)
                    count_2 = count_2 +1;
                    signal_confirmation(i,j) = 1;
                    D(i,j)= sqrt((xx(i)-xx(j))^2+(yy(i)-yy(j))^2);
                    D(j,i)= sqrt((xx(j)-xx(i))^2+(yy(j)-yy(i))^2);
                    ASHISH_MC(i,:) = 0;
                    STP_Matrix(j,:) = STP_Matrix(j,i)*STP_Matrix(i,:);
                    ASHISH_MC(j,:) = STP_Matrix(j,:);
                    ASHISH_MC(j,j) = 0;
                    for z = 1:number_of_nodes
                        if sum(signal_confirmation(:,z)) == 1
                            ASHISH_MC(:,z) = 0;
                       end
                    end
                   ASHISH_MC_2 = ASHISH_MC;
                   blah(count_2) = i;
                end
            end
        end
    end
    max_1 = max(max(ASHISH_MC));
end



H = graph(D);
H.Edges.Weight;

p = plot(H,'EdgeLabel',H.Edges.Weight);
p.XData = xx;
p.YData = yy;
for i=1:pointsNumber   
    for j=1:pointsNumber
     if R(i,j)==0
         highlight(p,[i,j],'EdgeColor', 'r')
     end
    end
end
