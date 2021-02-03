%Wing Properties 
d1=0.393;%Aerofoil length
d2=0.1; %Length to Lift Force
d3=0.165; %length to Centre of Mass
b=1.25; %aerofoil depth (NOT USED)
a=1.058; %Surface Area
m=4.830; %Mass of Aerofoil
theta=10; %Angle of Attack 
Ig=0.6;  %Moment of Inertia 
 
%Dynamic and Atmospheric Properties
Vi=26.82; %Car Velocity
Cl=10*theta; %Lift Coefficient
rho=1.225; %Density of air at sea level
 
%mesh grid Setup 
xx=0:10:15000;
yy=0:10:15000;
[k1,k2]=meshgrid(xx,yy);
 
Non_Oscillatory_k1=zeros(1500);
Non_Oscillatory_k2=zeros(1500);
stable_k1=zeros(1500);
stable_k2=zeros(1500);
flutter_k1=zeros(1500);
flutter_k2=zeros(1500);
 
%stability for-loop
for i=1:length(k1);
    for j=1:length(k2)
   k1_i=k1(i,j);
   k2_i=k2(i,j);
 
%State Space Model
a11=(k1_i+k2_i)/m;
a12=(((k1_i*d3)-(k2_i*(d1-d3)))/m)-((5*rho*a*(Vi^2))/m);
a21=(((k1_i*d3)-(k2_i*(d1-d3)))/Ig);
a22=(((k1_i*d3^2)+(k2_i*(d1-d3)^2))/Ig)-(((5*rho*a*(Vi^2))*(d3-d2))/Ig);

Amatrix=[a11 a12;a21 a22]; 
Mmatrix=[m 0;0 Ig];%Mass Matrix
Kmatrix=[(k1_i+k2_i) (k1_i*(d1-d3))-(k2_i-(d1-d3));(k1_i*d3)-(k2_i*(d1-d3)) ((k1_i*d3)^2)+((k2_i*(d1-d3))^2)];%Stiffness matrix
zero=zeros(2); %zero matrix
unit=eye(2); %unit matrix
Stability_Matrix=[zero unit;Amatrix zero]; %State Space matrix
Eigen_Values=eig(Stability_Matrix);%eigen values of the X_dot matrix
 
%Stability
A=((a11+a22 ))/2;
B=(a11*a22)-(a12*a21);
AA=A^2;
if (B<=0)
    disp('Non Oscillatory Unstable Divergence')
    Non_Oscillatory_k1(i,j)=k1_i;
    Non_Oscillatory_k2(i,j)=k2_i;
 
elseif (0<B)&&(B<=AA)
    disp('Stable Oscillation')
    stable_k1(i,j)=k1_i;
    stable_k2(i,j)=k2_i;
 
elseif (B>(AA))
    disp('Oscillatory Divergence "Flutter"')
 flutter_k1(i,j)=k1_i;
 flutter_k2(i,j)=k2_i;
    
end
    end
end
%Stability condition counters
No_of_non_Oscillatory=nnz(Non_Oscillatory_k1)
No_of_stable=nnz(stable_k1)
No_of_flutter=nnz(flutter_k1)
 
%Plotting commands 
hold on
plot(stable_k1,stable_k2,'green')
plot(Non_Oscillatory_k1,Non_Oscillatory_k2,'yellow')
plot(flutter_k1,flutter_k2,'red')
hold off
