clear
close all
clc

%%

%These are the reference results for the original point and the intial directions with vectors representing X,Y,Z axes
%of our intial frame when the robot is static without any translational or rotational inputs  
  %Point_Origin = [-370.673201266023,-304.185539730741,718.259481452566]; %Point of Origin for our frame
  %Q_Origin= [-1802.44593990606,223438.810653097,93696.9696140344];  % Q vector resmbles the X axis Direction
  %L_Origin= [-230845.344179285,4956.49219622259,-117033.595296806]; % L vector resmbles the Y axis Direction
  %S_Origin= [47315.0625741695,-72002.9479615230,-6075.56181031615]; % S vector resmbles the Z axis Direction

%Function to get points from Excel sheets (the number is based on the sheet number for each measured position)
%Number 6 for example means that we take the results from sheet 6 of our readings and so on
[num,txt,raw]= xlsread('GoniometerTest.xlsx',10);

%This is used to be an array for multiple positions but We make a simple change to become easier for the user using only one position each time with a predefined intial position
MeasuredPositions = 1;

%These are simple tables for gathering the input
P1=zeros(12,3);P2=zeros(12,3);P3=zeros(12,3);
Q1=zeros(12,3);Q2=zeros(12,3);Q3=zeros(12,3);
R1=zeros(12,3);R2=zeros(12,3);R3=zeros(12,3);

for i=1:size(MeasuredPositions,1)
    
    % XY-Plane function
    j=1:3;
    
    P1(i,:)=num(MeasuredPositions(i),j);
    P2(i,:)=num(MeasuredPositions(i)+1,j);
    P3(i,:)=num(MeasuredPositions(i)+2,j);
    XY_normal(i,:) = cross(P2(i,:)-P1(i,:), P3(i,:)-P1(i,:));
    
    d1(i,:) = dot(XY_normal(i,:),P1(i,:)); % d=-(ax0+by0+cz0)
    
    % XZ-Plane
    j=5:7;
    
    Q1(i,:)=num(MeasuredPositions(i),j);
    Q2(i,:)=num(MeasuredPositions(i)+1,j);
    Q3(i,:)=num(MeasuredPositions(i)+2,j);
    XZ_normal(i,:) = cross(Q2(i,:)-Q1(i,:), Q3(i,:)-Q1(i,:));
    
    d2(i,:) = dot(XZ_normal(i,:),Q1(i,:));  % d=-(ax0+by0+cz0)
    
    % YZ-Plane
    j=9:11;
    
    R1(i,:)=num(MeasuredPositions(i),j);
    R2(i,:)=num(MeasuredPositions(i)+1,j);
    R3(i,:)=num(MeasuredPositions(i)+2,j);
    YZ_normal(i,:) = cross(R2(i,:)-R1(i,:), R3(i,:)-R1(i,:));
    
    d3(i,:) = dot(YZ_normal(i,:), R1(i,:));   % d=-(ax0+by0+cz0)    
    
 
    A(i,:,:)=[YZ_normal(i,:);XY_normal(i,:);XZ_normal(i,:)];
    D(i,:)=[d3(i,:);d1(i,:);d2(i,:)];

    
    %Intersection of three planes is a PointOrigin
    
    PointOrigin(i,:)=squeeze(A(i,:,:))\squeeze(D(i,:)')   
      
    %X,Y,Z points before getting the othrogonal vectors 
    
    Z(i,:)=YZ_normal(i,:)- PointOrigin(i,:);
    Y(i,:)=XZ_normal(i,:)- PointOrigin(i,:);
    X(i,:)=XY_normal(i,:)- PointOrigin(i,:);
   
    %Develop Orthogonal vectors for acurate frames
    
    Q(i,:)=cross(X(i,:),PointOrigin(i,:)); %Orthogonal for X
    L(i,:)=cross(Y(i,:),PointOrigin(i,:)); %Orthogonal for Y
    S(i,:)=cross(Z(i,:),PointOrigin(i,:)); %Orthogonal for Z

 
   %Ploting Frame 
   
   
   PaintVector(Q(i,:),PointOrigin(i,:),'g') %X-Axis
   hold on
   PaintVector(L(i,:),PointOrigin(i,:),'b') %Y-Axis
   hold on
   PaintVector(S(i,:),PointOrigin(i,:),'r') %Z-Axis
   hold on 
    
end

%%
% The word previous after each parameter or variable points to the previous frame
[num,txt,raw]= xlsread('GoniometerTest.xlsx',1);

%This is used to be an array for multiple positions but We make a simple change to become easier for the user using only one position each time with a predefined intial position
MeasuredPositions = 1;

%These are simple tables for gathering the input
F1=zeros(12,3);F2=zeros(12,3);F3=zeros(12,3);
G1=zeros(12,3);G2=zeros(12,3);G3=zeros(12,3);
H1=zeros(12,3);H2=zeros(12,3);H3=zeros(12,3);

for o=1:size(MeasuredPositions,2)
    
    % XY-Plane function
    u=1:3;
    
    F1(o,:)=num(MeasuredPositions(o),u);
    F2(o,:)=num(MeasuredPositions(o)+1,u);
    F3(o,:)=num(MeasuredPositions(o)+2,u);
    XY_normal_previous(i,:) = cross(F2(o,:)-F1(o,:), F3(o,:)-F1(o,:));
    
    d1_previous(i,:) = dot(XY_normal_previous(o,:),F1(o,:)); % d=-(ax0+by0+cz0)
    
    % XZ-Plane
    u=5:7;
    
    G1(o,:)=num(MeasuredPositions(o),u);
    G2(o,:)=num(MeasuredPositions(o)+1,u);
    G3(o,:)=num(MeasuredPositions(o)+2,u);
    XZ_normal_previous(o,:) = cross(G2(o,:)-G1(o,:), G3(o,:)-G1(o,:));
    
    d2_previous(i,:) = dot(XZ_normal_previous(o,:),G1(o,:));  % d=-(ax0+by0+cz0)
    
    % YZ-Plane
    u=9:11;
    
    H1(o,:)=num(MeasuredPositions(o),u);
    H2(o,:)=num(MeasuredPositions(o)+1,u);
    H3(o,:)=num(MeasuredPositions(o)+2,u);
    YZ_normal_previous(o,:) = cross(H2(o,:)-H1(o,:), H3(o,:)-H1(o,:));
    
    d3_previous(o,:) = dot(YZ_normal_previous(o,:), H1(o,:));   % d=-(ax0+by0+cz0)    
    
 
    A_previous(o,:,:)=[XY_normal_previous(o,:);XZ_normal_previous(o,:);YZ_normal_previous(o,:)];
    D_previous(o,:)=[d1_previous(o,:);d2_previous(o,:);d3_previous(o,:)];

    
    % Intersection of three planes is a PointOrigin
    
    PointOrigin_previous(o,:)=squeeze(A_previous(o,:,:))\squeeze(D_previous(o,:)')   
      
    %X,Y,Z points before getting the othrogonal vectors 
    
    Z_previous(o,:)=YZ_normal_previous(o,:)- PointOrigin_previous(o,:);
    Y_previous(o,:)=XZ_normal_previous(o,:)- PointOrigin_previous(o,:);
    X_previous(o,:)=XY_normal_previous(o,:)- PointOrigin_previous(o,:);
   
    %Develop Orthogonal vectors for acurate frames
    
    Q_previous(o,:)=cross(X_previous(o,:),PointOrigin_previous(o,:)); %Orthogonal for X
    L_previous(o,:)=cross(Y_previous(o,:),PointOrigin_previous(o,:)); %Orthogonal for Y
    S_previous(o,:)=cross(Z_previous(o,:),PointOrigin_previous(o,:)); %Orthogonal for Z
    
   %Ploting Frame    
   PaintVector(Q_previous(o,:),PointOrigin_previous(o,:),'m') %X-Axis
   hold on
   PaintVector(L_previous(o,:),PointOrigin_previous(o,:),'y') %Y-Axis
   hold on
   PaintVector(S_previous(o,:),PointOrigin_previous(o,:),'k') %Z-Axis
   hold on 
end
%%
%Rotation Vectors to calculate the 3x3 rotation matrix
    Alpha = vrrotvec(Q(i,:),Q_previous(o,:)) %Rotation between X axes of the previous frame to the next frame
    Beta  = vrrotvec(L(i,:),L_previous(o,:)) %Rotation between Y axes of the previous frame to the next frame
    Gamma = vrrotvec(S(i,:),S_previous(o,:)) %Rotation between Z axes of the previous frame to the next frame
%%

%DH Table constants and variables

%Variables based on the results from the optics lab
d2=0;
d3=0;
d4=0;
Theta5=0;
Theta6=0;
Theta7=0;
%Constant lenghts in mm from Inventor 
L1=624.643;
L2=255.5;
L3=101.01;
L4=129.307;
L5=109.500;
L6=52.994;
L7=41.88;
%Forward Kinematics Matrix
Forward_Kinematics_Matrix = [-cos(Theta6)* sin(Theta7)	-cos(Theta6)*cos(Theta7)	sin(Theta6) 	d4+L1+L3+(L5*sin(Theta6))+(L7*sin(Theta6))-(L6*cos(Theta6)*sin(Theta7));
(cos(Theta5)*cos(Theta7))-(sin(Theta5)*sin(Theta6)*sin(Theta7))	 -(cos(Theta7)*sin(Theta5)*sin(Theta6))-(cos(Theta5)*sin(Theta7))	-cos(Theta6)*sin(Theta5)	d2+L4*sin(Theta5)-L5*cos(Theta6)*sin(Theta5)-L7*cos(Theta6)*sin(Theta5)+L6*cos(Theta5)*cos(Theta7)-sin(Theta5)*sin(Theta6)*sin(Theta7);
(cos(Theta7)*sin(Theta5))+(cos(Theta5)*sin(Theta6)*sin(Theta7)) 	(cos(Theta5)*cos(Theta7)*sin(Theta6))-(sin(Theta5)*sin(Theta7))	 cos(Theta5)*cos(Theta6)	d3+L2-L4*cos(Theta5)+L5*cos(Theta5)*cos(Theta6)+L7*cos(Theta5)*cos(Theta6)+L6*cos(Theta7)*sin(Theta5)+cos(Theta5)*sin(Theta6)*sin(Theta7);
0	0	0	1]

%X,Y and Z from forward kinematics transformation matrix
Xf= d4+L1+L3+L5*sin(Theta6)+L7*sin(Theta6)-L6*cos(Theta6)*sin(Theta7);
Yf= d2+L4*sin(Theta5)-L5*cos(Theta6)*sin(Theta5)-L7*cos(Theta6)*sin(Theta5)+L6*cos(Theta5)*cos(Theta7)-sin(Theta5)*sin(Theta6)*sin(Theta7);
Zf= d3+L2-L4*cos(Theta5)+L5*cos(Theta5)*cos(Theta6)+L7*cos(Theta5)*cos(Theta6)+L6*cos(Theta7)*sin(Theta5)+cos(Theta5)*sin(Theta6)*sin(Theta7);

%%
Rotation_Matrix = [Alpha; Beta; Gamma; 0 0 0 1] ; %Rotation Matrix is out based on the rotation between the different frames
PointOrigin_Transpose = [transpose(PointOrigin);0]; %The transpose is occured to replace the unneeded coloumn of the rotation matrix
Extract_Not_Needed_Coloumn = Rotation_Matrix(:, 4); %The needless coloumn is out
Correct_Coloumn=Extract_Not_Needed_Coloumn+PointOrigin_Transpose; 
Rotation_Matrix(:,4) = Correct_Coloumn; %The correct coloumn is introduced to the rotation matrix
Inverse_Kinematics_Matrix= Rotation_Matrix %The inveres kinematics matrix is produced after adding the correct coloumn to the rotation matrix
Error_Matrix= Inverse_Kinematics_Matrix-Forward_Kinematics_Matrix %Error matrix is defined as the difference between the forward and inverse kinematics matrices for the same frame or position
