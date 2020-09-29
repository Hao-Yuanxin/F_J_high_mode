clear;
%   This program calculates the dispersion curve of Love wave
%   reference： Xiaofei Chen(1993)，A systematic and efficient method of 
%   computing normal modes for multilayered half-space
%   Geophys. J. Int. (1993) 115,391-409
%   structure  parameters
global N z_bottom rho beta alpha miu Epsilon;
global delta_c delta_samll_c c_love_min c_love_max  f_all mode_number;      % Search step size ,cmax, cmin

% size : z_bottom (N) 
% size : z_bottom rho beta alpha miu (N+1) 
z_bottom=[18.0,24.0,30.0];  % depth to bottom,km,last layer is infinite
rho=[2.8,2.9,3.1,3.3];          % density for each layer,g/cm3
beta=[3.50,3.65,3.90,4.70];     % velocity for S wave,km/s
alpha=[6.00,6.30,6.70,8.20];    % velocity for P wave,km/s
N=length(z_bottom);         % the number of layers,half-space is N+1;
miu=beta.^2.*rho;               % u，Lame coefficient 
Epsilon = 1*10^(-6);   %  A small number， control error precision
delta_c = 0.001;  % Search step size of c
delta_samll_c = 0.0001;
       
c_love_min = 3.0;      % Set the minimum velocity for love wave
c_love_max = 5.5;      % Set the minimum velocity for love wave
mode_number = 40;      % number of mode
f_all = 0.1:0.1:5;       % The range of frequency



%SH case  
% meanings of the parameters, see Xiaofei Chen（1993）
c=c_love_min:delta_c:c_love_max;      % love wave velocity


%%  T = 1s      plot secular function ,c is  
% plot Love wave's secular function at T = 1s 
w=2*pi*1;          %T =1s w=2*pi*f
Nc = length(c);      %Nc length of c[];
V_secu = zeros(1,Nc);  % secular function f(c)
V_secu_R = zeros(1,Nc); % the real part of V_secu
V_secu_I = zeros(1,Nc); % the imaginary part of V_secu


for m=1:Nc
    V_secu(m)   = SECULAR(w,c(m));
    V_secu_R(m) = real(V_secu(m));
    V_secu_I(m) = imag(V_secu(m));
end
figure(1)
plot(c,V_secu_R,'k-',c,V_secu_I,'r-.');

xlabel('Velocity(km/s)');
ylabel('Secular function');
title('Secular function T = 1s');
hold on
root_all_1s = findroot_all(w);
N_1s = length(root_all_1s);
y_1s = zeros(1,N_1s);
plot(root_all_1s,y_1s,'o','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b');
legend('Real','Imaginary','root');
hold off

%%
 



figure(2)
file_name ='~/Subject/CBV_highmode/Forward/dispersion_curves/Love_dispersion_Chen.txt';

plot_write_Love_dispersion(file_name);



                                                                                





%%
function plot_write_Love_dispersion(file_name)
    global N z_bottom rho beta alpha ;
    global f_all c_love_min c_love_max mode_number ;
  
    f = fopen(file_name,'w+');
    
    % write information
    fprintf(f,'!\t \t Love wave''s dispersion curves \n' );
    fprintf(f,'!\t \t   Parameter setting \n' );
    fprintf(f,'! There are %d + 1 interfaces in total;the free interface is the zero interface\n',N);
    fprintf(f,'! depth of each interface(except 0 interface) ');
    for m = 1: N
        fprintf(f,'%.2f ',z_bottom(m));
    end
    fprintf(f,'\n! Density of each layer                       ');
    for m = 1: N+1
        fprintf(f,'%.2f ',rho(m));
    end
    fprintf(f,'\n! alpha（P wave velocity) of each layer       ');
    for m = 1: N+1
        fprintf(f,'%.2f ',alpha(m));
    end
    fprintf(f,'\n! beta （S wave velocity) of each layer       ');
    for m = 1: N+1
        fprintf(f,'%.2f ',beta(m));
    end
    N_f = length(f_all);
    fprintf(f,'\n! frequency     range : f_min = %.2f, f_max = %.2f ', f_all(1),f_all(N_f));
    fprintf(f,'\n! Love velocity range : c_min = %.2f, c_max = %.2f ', c_love_min, c_love_max );
    fprintf(f,'\n! mode number : %d (0 to %d)', mode_number ,mode_number-1 );
    
    
    % polt and write Love dispersion
    fprintf(f,'\n! The first row is frequency（Hz) ');
    fprintf(f,'\n! The second row is love wave''s velocity \n ');
    w_all = f_all*2*pi;
    V_M = get_V_M(w_all);
    [N_f,N_mode] = size(V_M);    % N_w,number of w;N_mode,number of mode
    for m = 1:N_mode   %each column, each mode
         fprintf(f,'#  \t %d \n',m-1);
         f_temp = [];
         c_temp = [];
         for n = 1:N_f  % one mode  
             if (V_M(n,m) ~= 0)     % Find the starting point of dispersion curve
                 f_temp = [f_temp,f_all(n)];
                 c_temp = [c_temp,V_M(n,m)];
             end
             
         end
         fprintf(f,'%.4f ',f_temp); % write frequency row
         fprintf(f,'\n');
         fprintf(f,'%.4f ',c_temp);
         fprintf(f,'\n');
         plot(f_temp,c_temp,'k-')
         hold on
         
    end
                 

   

end
























function  lambda_d = Lambda_d(z,zj_1,v_j)     % j=1,2...N (14a)

    lambda_d = exp(-v_j*(z-zj_1));

end

function  lambda_u = Lambda_u(z,zj,v_j)       % j=1,2...N (14a)

    lambda_u = exp(-v_j*(zj-z));

end

function secular_w_c = SECULAR(w,c)     
    % secular function vlaue at (w,c) 
    global N z_bottom  beta miu  Epsilon;
    % secular function at single frequency and single velocity
    E=zeros(2,2,N+1);      %3D array,3th dimension represents the number of layers
    nu=zeros(1,N);       %  v,    (k^2-(w/beta)^2)^(1/2)

    RT_matrix=zeros(2,2,N-1); % Reflection transmission matrix  （17a)
    RT_N=zeros(2,1);          % (Td Rdu)T (17b)
    Td  = zeros(1,N);            
    Rud = zeros(1,N);
    Rdu = zeros(1,N);
    Tu  = zeros(1,N);
    Rdu_G = zeros(1,N);  % generalized reflection and transmission coefficients
    Td_G  = zeros(1,N); 
 

% when c == beta(j) ,the secular function is singular .in this case, i plus
% a small number to c
    for n= 1:N+1
        if c == beta(n)
            c = c+Epsilon;
        end
    end
    for m = 1: N+1
        k=w/c;
        nu(m)=(k^2-(w/beta(m))^2)^(1/2);
        % E for SH 
        E(1,1,m)=1;                 %   (11)
        E(1,2,m)=1;                 %   (11)
        E(2,1,m)=-miu(m)*nu(m);     %   (11)
        E(2,2,m)=miu(m)*nu(m);      %   (11)
    end
    clear m


    for m = 1:N-1
        if m == 1                   %  z_0==0;
            lmd = Lambda_d(z_bottom(m),0,nu(m));     % lambda11 (17a)
        else 
            lmd = Lambda_d(z_bottom(m),z_bottom(m-1),nu(m)); % lambda11 (17a)
        end
        lmu = Lambda_u(z_bottom(m),z_bottom(m+1),nu(m+1));   % lambda22 (17a)
        RT_matrix(:,:,m)=[E(1,1,m+1),-E(1,2,m);E(2,1,m+1),-E(2,2,m)]\ ...
        [E(1,1,m),-E(1,2,m+1);E(2,1,m),-E(2,2,m+1)]*[lmd,0;0,lmu] ;  % (17a)  
    %   inv(A)*B == A\B
        Td(m)  = RT_matrix(1,1,m);
        Rud(m) = RT_matrix(1,2,m);
        Rdu(m) = RT_matrix(2,1,m);
        Tu(m)  = RT_matrix(2,2,m);
    
    end

% j= N, interface
    lmd_N = Lambda_d(z_bottom(N),z_bottom(N-1),nu(N));
    RT_N  = [E(1,1,N+1),-E(1,2,N);E(2,1,N+1),-E(2,2,N)]\[E(1,1,N)*lmd_N;E(2,1,N)*lmd_N];
%   inv(A)*B == A\B
    Td(N) = RT_N(1,1);
    Rdu(N)= RT_N(2,1);

% generalized reflection and transmission coefficients (Rdu_G Td_G)
% j=N, interface
    Rdu_G(N) = Rdu(N);
    Td_G(N)  =Td(N);
% j=1,2,3...N-1 interface
    for m = N-1:-1:1
        Td_G(m) = (1-Rud(m)*Rdu_G(m+1))\Td(m); % （21） note! it is \ ,not / !!!
        Rdu_G(m)= Rdu(m)+Tu(m)*Rdu_G(m+1)*Td_G(m);
    end

% Rud_G(0)
    Rud_G0  = -E(2,1,1)\E(2,2,1)* Lambda_u(0.0,z_bottom(1),nu(1));
    secular_w_c = 1-Rud_G0*Rdu_G(m);


end




function [root_one,x_right]=findroot_1(w0,c_start,step) 
% x_left is the right boundary of root interval
% case : w = w0,find a sinale root
    global c_love_max  Epsilon;
    
    %  Find the interval of root
    F0 = imag(SECULAR(w0,c_start));
    c_temp = c_start + step;
    F1 = imag(SECULAR(w0,c_temp)); 
    
    while (F0/F1 > 0)
       
        if (c_temp > c_love_max)
            root_one = [];
            x_right = [];
            return
        else
            c_temp = c_temp +step;
            F0 = F1;
            F1 = imag(SECULAR(w0,c_temp));

        end
    end
    x_right = c_temp;
    % c_temp must less than B 
    if (c_temp >= c_love_max) % must > = ,if not ,may cause x1 = x2，the answer is Inf
        c_temp = c_love_max;
        F1 = imag(SECULAR(w0,c_temp));
        if (F0/F1 > 0)  %Judge whether the interval has a solution
            root_one = [];
            x_right = [];
            return
        end    
        
    end
   
    % The interval of the root is found
    
    % find the root using secant method。
    m = 0;      
    n_iter = 500;   %Number of iterations
    x1 = c_temp-step;
    x2 = c_temp;
    while (m<n_iter)
        m = m+1;
        x_temp = x1-F0/(F1-F0)*(x2-x1);
      
        if (abs(x_temp-x2) > Epsilon)
            x1 = x2;
            x2 = x_temp;
            F0 = imag(SECULAR(w0,x1)); 
            F1 = imag(SECULAR(w0,x2)); 
        else
            break
        end
    end
    root_one =  x_temp;
end

function c_begin = find_begin(w0)
    global    beta delta_samll_c ;
    
    N = length(beta);
    c1 = beta(1);
   
    for m = 2:N
        c2 = beta(m);
        if (c1 > c2)
            c1 = c2;
        end
    end
    c_temp = c1;
  
    step = delta_samll_c;
    while(1)      % Infinite cycle    
        F = imag(SECULAR(w0,c_temp)); 
        if abs(F) > 0.1
            break
        end
        c_temp =  c_temp+step;
    end
    c_begin = c_temp; 
end

function root_all = findroot_all(w0)
    % find all roots when w = w0
    global c_love_max  delta_c delta_samll_c Epsilon;
    root_all = [];
    c_start = find_begin(w0);
    
    % first root
    step = delta_samll_c;
    %{
    [c_temp ,c_right] = findroot_1(w0,c_start,step); % root of imaginary
    if (isempty(c_right))
         return
             
    end  
    if (abs(c_temp-c_right)< Epsilon *10)
        c_right = c_right + Epsilon*10;
    end
    if (c_right>=c_love_max)
         if real(SECULAR(w0,c_temp)) < 0.3    % judge real part
             root_all = [root_all,c_temp];
         end

         return
    end

    
    
    if real(SECULAR(w0,c_temp)) < 0.3    % judge real part
             root_all = [root_all,c_temp];
  
    end
    %}
    m = 0;
     while(m <= 5 )
         m = m+1;
         [c_temp ,c_right] = findroot_1(w0,c_start,step); % root of imaginary  
         if (isempty(c_right))
            return
             
         end
         
         if (abs(c_temp-c_right)<Epsilon*10)      % case: c_temp == c_right
             c_right = c_right + Epsilon*10;
         end
         if (c_right>=c_love_max)
             if real(SECULAR(w0,c_temp)) < 0.3    % judge real part
                 root_all = [root_all,c_temp];  
             end
             
             return
         end
         c_start = c_right;
         if real(SECULAR(w0,c_temp)) < 0.3    % judge real part
              root_all = [root_all,c_temp];
  
         end
    end
    c_start = c_right;
    
    step = delta_c;
    
    
    while(1)
         
         [c_temp ,c_right] = findroot_1(w0,c_start,step); % root of imaginary  
         if (isempty(c_right))
            return
             
         end
         
         if (abs(c_temp-c_right)<Epsilon*10)      % case: c_temp == c_right
             c_right = c_right + Epsilon*10;
         end
         if (c_right>=c_love_max)
             if real(SECULAR(w0,c_temp)) < 0.3    % judge real part
                 root_all = [root_all,c_temp];  
             end
             
             return
         end
         c_start = c_right;
         if real(SECULAR(w0,c_temp)) < 0.3    % judge real part
              root_all = [root_all,c_temp];
  
         end
    end
end
function V_M = get_V_M(w_all)
% w is an array of angular frequencies
    global mode_number;
    N_w = length(w_all);   % length of w_all
    V_M = zeros(N_w, mode_number);  % Velocity matrix，Each row corresponds to a frequency,
    %   each column corresponds to a mode
    
    for m = 1:N_w
        root_temp = findroot_all(w_all(m)); 
        N_root = length(root_temp);
        if (N_root > mode_number)     % The number of columns cannot exceed mode_number
             N_root = mode_number;    
        end
        for n = 1:N_root
            V_M(m,n) = root_temp(n);
        end

    end

    
    
    
    
end






