
function Dicplacment =truss_4

%........ Assumptions ........
%- E = constant .
%- A = constant for all members
%- input supports after free nodes (arange nodes according its degree of freedom )
%- joint can have any number of Load 
%- angle (Degree )// other parameters uniteless
%- no settelment consideration



%input Script (Preprocess)----------------------------

No_members=input('Number of members =');

%Modelas of Elastisty input 
% con used as condation
con=input('IF E (Modelas of Elastisty) is constant press (0) , variabel (1)');
 E=zeros(No_members,1);
if con == 1
   
    for e=1:No_members
    fprintf('For member   %d',  e  )
    E(e,1)=input('Modelas of Elastisty =');
    end
end
    
if con==0 
   E(1:end,1)=input('Modelas of Elastisty = ');  
end
%end of E input

%input Script (Preprocess)----------------------------
                                         
Area=input('corss section Area = ');

No_nodes=input('Number of nodes =');
%-----------------------------------------------
x=zeros(1,No_nodes); %creataion of cordinates Matrix (1 Row Matrix)
y=zeros(1,No_nodes); 
for e=1:No_nodes %input the cordinates of each node (for loop)
                 % {e} is the counter  
   fprintf('input cordinates for Node : %d ',e)
   x_1=input ('  X  =');
   x(:,e)=x_1+ x(:,e);
   fprintf('input cordinates for Node : %d ',e)
   y_1=input (' y  =' );
   y(:,e)=y_1+ y(:,e);
end  %end of input the cordinates of Nodes
%.......................
  %type of Nodes
  ('input type of each node')
  R=zeros(1,No_nodes); %Rigidty Matrix for each node (1 Row Matrix)
  %R=[0 1 2 2 2 0]
  (' (Hinge support = 0 , Roller support = 1 , Free = 2) ')
for e=1:No_nodes
   fprintf('for Node: %d ',e)
    R1=input('=');
   R(:,e)=R1+R(:,e);
end
%.........................
 %matching members 
 J=zeros(No_members,2);
 % J = [start   End]
 %      start   End
 ('Start & End Joint for each member')
 for e=1:No_members
     fprintf('Member %d  ',e);
     J1=input('  start : ');
     J(e,1)=J(e,1)+J1;
      J1=input('End : ');
     J(e,2)=J(e,2)+J1;
 end
 %........................................
 %input Loads..............
 %can put any number of loads on Joint
 fprintf('Input Number of forces')
 No_loads=input('=');
 loads=zeros(No_loads,3);
 %.. loads=[ No.Joint   Value of load    angle ]
 ('[ Joint Nubmer   Value of load    angle of load ]')
 %%%%
 for e=1:No_loads
    No.Joint=input('Joint Nubmer  = ');
    loads(e,1)=loads(e,1)+No.Joint;
    %............
    value=input('Value of load = ');
    loads(e,2)=loads(e,2)+value;
    %............
    angel=input('angle of load = ');
    loads(e,3)=loads(e,3)+angel;
 end
 loads;
%.............End of preprocess.........................
%........................................................................................
% frist itration of solveing......................
%lengh and angle for each member ..

% Member matrix [length   inclination]............
members=zeros(No_members,2);
for e=1:No_members
    A=J(e,1);
    B=J(e,2);
    x1=x(1,A);
    x2=x(1,B);
    y1=y(1,A);
    y2=y(1,B);
    dx=x2-x1; 
    dy=y2-y1;
   members(e,1)=sqrt(dx^2 + dy^2 ) ; % Length 
   members(e,2)=atand(dy/dx);% angle = tan(Diff Y / Diff X)
end 
('[length   inclination]')
members
%............................
% second Itration of solving.............
%  K local ....
Ks=zeros(No_nodes*2,No_nodes*2);% K structure 
% Kl=[1 0 -1 0
%     0 0 0 0
%     -1 0 1 0
%     0 0 0 0];
 
for e=1:No_members
   % k=Kl*(E*A/members(e,1));
    %taransformation.............
    c=cosd(members(e,2));
    s=sind(members(e,2));
    cs=s*c;
    c2=c*c;
    s2=s*s;
    Transformation=[c2 cs -c2 -cs
                    cs s2 -cs -s2   
                   -c2 -cs c2 cs
                   -cs -s2 cs s2]
               L2=members(e,1)
               kk =(E(e,1)*Area)/ L2 
             fprintf('%d',e)%No. of member
    Kg =  (Transformation * ( kk )) % K global for element
    A=J(e,1);
    B=J(e,2);         % use A*2-1 & A*2 to access the right cell in Ks becase the count of matrix cells start with 1
    
    Ks(A*2-1:A*2,A*2-1:A*2)=Ks(A*2-1:A*2,A*2-1:A*2)+Kg(1:2,1:2);%Diagonal partations
    Ks(B*2-1:B*2,B*2-1:B*2)=Ks(B*2-1:B*2,B*2-1:B*2)+Kg(3:4,3:4);
    
    Ks(A*2-1:A*2,B*2-1:B*2)=Ks(A*2-1:A*2,B*2-1:B*2)+Kg(1:2,3:4);%symetric partations
    Ks(B*2-1:B*2,A*2-1:A*2)=Ks(B*2-1:B*2,A*2-1:A*2)+Kg(3:4,1:2);
    
    % submartises concept is used^_^ its the  most powerful thing in Matlab
   
end  
Ks
%..... Third Itration of solving.............
%Dividing KS [ Kuu   Kur]........
%              Kru   Krr


%we need to count the dimantion of sub matrices 
%use For loop to count the degrees of freedom
%use Q to Kuu --- w to Kur --- r to Kru --- t to Krr in the fact we didn't need to calculate any thing except Q ... the others well be the rest of the matrix
%------------to consider settlement rescript this part
Q=0;
for e=1:No_nodes
    if R(1,e)==2
      Q=Q+2;  
    elseif R(1,e)==1
      Q=Q+1;
    end 
end
Kuu=Ks(1:Q,1:Q) %use Kuu to calculate Displacment
Kru=Ks(Q+1:end,1:Q)%use Kru to calculate Reactions 

F1=zeros(Q,1);%force matrix .. one column
% to get Number of loads 
n=size(loads);
% use For loop to inputing loads in force matrix
% load in any direction
for e=1:n(1,1)
% if loads(e,3)==0
     F1((loads(e,1)*2)-1,1)=F1((loads(e,1)*2)-1,1)+(loads(e,2)*cosd(loads(e,3)));
 %elseif loads(e,3)==90
     F1((loads(e,1)*2),1)=F1((loads(e,1)*2),1)+(loads(e,2)*sind(loads(e,3)));
% end
end

F1
d1=Kuu^-1 *F1  % delta (Displacment) according to Dgree of fredom of every Node
F2=Kru*d1

%-----------------------------------------------------------------------------
%calculating internal forces 
%[FL]=[KL]*[T]t*[Dg]
%"we will make a big loop to calculate each variable and internal force for each member in the same ittraion"
%
for e=1:No_members
   
    %taransformation.............
    c=cosd(members(e,2));
    s=sind(members(e,2));
    
    Tt=[ c  s  0  0
        -s  c  0  0   
         0  0  c  s
         0  0 -s  c];
   %------------------------- 
   % K Local ..........
               L2=members(e,1);
               KL =((E(e,1)*Area)/ L2)*[ 1 0 -1 0 
                                    0 0  0 0
                                   -1 0  1 0
                                    0 0  0 0];
   %---------------------------
   % D global for each member 
   Dg=zeros(4,1);
    A=J(e,1);
    B=J(e,2);
    Ld1=length(d1);
  %  for i=1:4
        if (A*2)-1 <= Ld1
          Dg(1,1)=Dg(1,1)+ d1((A*2)-1 ,1);
        end    
            %..........
        if (A*2) <= Ld1
          Dg(2,1)=Dg(2,1)+ d1((A*2) ,1);
        end    
            %..........
        if (B*2)-1 <= Ld1
          Dg(3,1)=Dg(3,1)+ d1((B*2)-1 ,1);
        end    
            %..........
        if (B*2) <= Ld1
          Dg(4,1)=Dg(4,1)+ d1((B*2) ,1);
        end    
            %..........
            
  %  end
      % end of calculating Dg ---------
      Dg
   
   %solve equation  [FL]=[KL]*[T]t*[Dg]
   fprintf('%d',e)%No. of member
   FL = KL * Tt * Dg
   
                
               
end