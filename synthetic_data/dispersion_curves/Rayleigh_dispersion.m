clear;

%   Author: 郝渊新 Hao Yuanxin 
%   Company: IGG
%   email: haoyuanxin19@mails.ucas.ac.cn
%   This program calculates the dispersion curve of Rayleigh wave
%   reference： Xiaofei Chen(1993)，A systematic and efficient method of 
%   computing normal modes for multilayered half-space
%   Geophys. J. Int. (1993) 115,391-409
%   structure  parameters
global N z_bottom rho beta alpha miu Epsilon beta_min;
global delta_c delta_samll_c c_Ray_min c_Ray_max  f_all mode_number cr;      % Search step size ,cmax, cmin

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
       
c_Ray_min = 3.0;      % Set the minimum velocity for love wave
c_Ray_max = 5.5;      % Set the minimum velocity for love wave
mode_number = 40;      % number of mode
T = 100:-1:11;
f1 = 1./T;
f2 = 0.1:0.1:5;
f_all = [f1,f2] ;    % frequency_all





beta_min = beta(1);    % get the beta_min
for m=1:N
    if (beta_min>beta(m+1))
        beta_min = beta(m+1);
    end
end
        







%%  T = 1s      plot secular function ,c is  
% plot Love wave's secular function at T = 1s 

cr = getCR();   % w --> infinite ,love wave velocity
if c_Ray_min > cr
    fprintf('c_Ray_min is larger than Cr, program is terminated \n');
    return
end
%fprintf('c_fund = %.8f  c_r = %.8f ', c_fund,cr) ;

c  =  c_Ray_min:delta_c:c_Ray_max;
w=2*pi/31;                   %T =1s w=2*pi*f
Nc = length(c);             %Nc length of c[];
V_secu = zeros(1,Nc);       % secular function f(c)
V_secu_R = zeros(1,Nc);     % the real part of V_secu
V_secu_I = zeros(1,Nc);     % the imaginary part of V_secu


for m=1:Nc
    if (c < beta_min)
         V_secu_R(m)  = SECULAR_Ray(w,c(m));
         V_secu_I(m)  = 0;
    else 
        V_secu(m)   = SECULAR_Ray(w,c(m));
        V_secu_R(m) = real(V_secu(m));
        V_secu_I(m) = imag(V_secu(m));
    end
end
figure(1)
plot(c,V_secu_R,'k-',c,V_secu_I,'r-.');

xlabel('Velocity(km/s)');
ylabel('Secular function');
title('Secular function T = 10s');
hold on
%{
root_fund = root_find_appr(w,3.3); 
plot(root_fund,0,'o');
legend('Real','Imaginary','root');

%}
root_real=findroot_real(w);     % ture root
if (root_real)
    plot(root_real,0,'o');
end
% legend('Real','Imaginary','rootfund','roothigh');

%}
%root_fund = root_find_appr(w);
%{
if (root_fund)       % approximate root
    plot(root_fund,0,'o','MarkerFaceColor','b');
end
%}
%c_fund = c_fundamental_appr(w);
%plot(c_fund,0,'o');
%fprintf('%.8f',c_fund);
%
root_all_1s = findroot_all(w);
N_1s = length(root_all_1s);
y_1s = zeros(1,N_1s);
plot(root_all_1s,y_1s,'o','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b');
%}
 
%{
root_temp_real = findroot_real(w);  % c < beta_min
root_temp_imag = findroot_all(w);   % beta_max < c < Ray_max

root_temp = [root_temp_real,root_temp_imag];
root_fund = root_find_appr(w); 
root_temp = [root_fund,root_temp];
y_ns = zeros(1,length(root_temp));
plot(root_temp,y_ns,'o');
hold off
%}
%%
% w_all = 2*pi*f_all;
% V_M = get_Ray_V_M(w_all) ; 


figure(2)
file_name ='~/Subject/CBV_highmode/Forward/dispersion_curves/Ray_dispersion_Chen.txt';



plot_write_Ray_dispersion(file_name);








%%



function plot_write_Ray_dispersion(file_name)
    global N z_bottom rho beta alpha ;
    global f_all c_Ray_min c_Ray_max  mode_number ;
  
    f = fopen(file_name,'w+');
    
    % write information
    fprintf(f,'!\t \t Rayleigh wave''s dispersion curves \n' );
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
    fprintf(f,'\n! frequency     range : f_min = %.4f, f_max = %.4f ', f_all(1),f_all(N_f));
    fprintf(f,'\n! Rayleigh velocity range : c_min = %.4f, c_max = %.4f ', c_Ray_min, c_Ray_max );
    fprintf(f,'\n! mode number : %d (0 to %d)', mode_number ,mode_number-1 );
    
    
    % polt and write Love dispersion
    fprintf(f,'\n! The first row is frequency（Hz) ');
    fprintf(f,'\n! The matrix is love wave''s velocity \n');
    fprintf(f,'#      frequency（Hz）  \n');
    fprintf(f,'%.4f ',f_all);
    fprintf(f,'  \n');
    fprintf(f,'# Rayleigh wave''s velocity，each column corresponds to a mode number \n');
  
    w_all = f_all*2*pi;
    V_M = get_Ray_V_M(w_all);
    
    [N_f,N_mode] = size(V_M);    % N_w,number of w;N_mode,number of mode
    for m = 1:N_f
         fprintf(f,'%.4f ',V_M(m,:));
         fprintf(f,'\n');
    end
    for m = 1:N_mode   %each column, each mode
     
         f_temp = [];
         c_temp = [];
         for n = 1:N_f  % one mode  
             if (~ isnan(V_M(n,m)))     % Find the starting point of dispersion curve
                 f_temp = [f_temp,f_all(n)];
                 c_temp = [c_temp,V_M(n,m)];
             end
             
         end
         
         plot(f_temp,c_temp,'k-')
         hold on
         
    end
                 
    ylabel('Velocity(km/s)');
    xlabel('frequency /Hz');
    title('Dispersion curves of Rayleigh wave');
   

end

function root_fund=root_find_appr(w,root_last)  
    % Find roots in very narrow interval,
    % The upper limit of this interval is  the  root of the last order
    
    % step =10^(-12);
    zero_14 = 10^(-8);    % zero
    n_step = 1000 ;         % number of step
    c_fund = c_fundamental_appr(w);  % approximate fundamental Rayleigh wave velocity
    c_min = c_fund;    % Search interval
    c_max = root_last;
    if (abs(c_fund-root_last)<zero_14)
        root_fund = c_fund;
        return
    end
        
    step = abs(c_max-c_min)/n_step ;
    c_max = root_last+step;
    c_temp  = c_min-step;
    F0 = real(SECULAR_Ray(w,c_temp));
    c_temp  =  c_temp +step;
    F1 = real(SECULAR_Ray(w,c_temp)); 
    
    while (~(F0<0 && F1>0))
       
        if (c_temp > (c_max+step))
            root_fund = [];
            fprintf('no root in  narrow Search interval,f is %.2f \n',w/2/pi);
            return
        else
            c_temp = c_temp +step;
            F0 = F1;
            F1 = real(SECULAR_Ray(w,c_temp));

        end
    end
   
  

    % c_temp must less than  c_love_max
    if (c_temp >= (c_max+step)) % must > = ,if not ,may cause x1 = x2，the answer is Inf
        c_temp =  c_max+step;
        F1 = real(SECULAR_Ray(w,c_temp));
        if (F0/F1 > 0)  %Judge whether the interval has a solution
            root_fund = [];
            fprintf('no root in  narrow Search interval,f is %.2f \n',w/2/pi);
            return
        end    
        
    end
   
    % The interval of the root is found
    
    % find the root using secant method。
    m = 0;      
    n_iter = 500;   %Number of iterations
    x_min = c_temp-step;
    x1 = x_min;
    x2 = c_temp;
    while (m<n_iter)
        m = m+1;
        x_temp = x1-F0/(F1-F0)*(x2-x1);
        if (x1<x_min)
            x1 = x_min;
        elseif (x2<x_min)
            x2 = x_min;
        end
            
        
        if (abs(x_temp-x2) >  zero_14)
            x1 = x2;
            x2 = x_temp;
            F0 = real(SECULAR_Ray(w,x1)); 
            F1 = real(SECULAR_Ray(w,x2)); 
       
            
        else
            break
        
        end
    end
    if (abs(imag(SECULAR_Ray(w,x_temp))<10^(-3)))
        root_fund =  x_temp;
    else
        root_fund = [];
    end
    
end
    
    
    
    

function V_M = get_Ray_V_M(w_all)  % get the velocity matrix
% w is an array of angular frequencies
    global mode_number;
   
    N_w = length(w_all);   % length of w_all
    V_M = zeros(N_w, mode_number)*NaN;  % Velocity matrix，Each row corresponds to a frequency,
    %   each column corresponds to a mode
    
    for m = 1:N_w
        root_temp_real = findroot_real(w_all(m));  % c < beta_min
        root_temp_imag = findroot_all(w_all(m)) ; % beta_max < c < Ray_max
        root_temp = [root_temp_real,root_temp_imag];
      
        if (m>1 && (V_M(m-1,1) < root_temp(1)))
            root_fund = root_find_appr(w_all(m),V_M(m-1,1)); 
            root_temp = [root_fund,root_temp];
        end
        N_root = length(root_temp);
        if (N_root > mode_number)     % The number of columns cannot exceed mode_number
             N_root = mode_number;    
        end
        for n = 1:N_root
            V_M(m,n) = root_temp(n);
        end

    end

    
    
    
    
end



function root_real=findroot_real(w)
    global  delta_c beta beta_min cr N  Epsilon ;
    root_real = [];
    c = cr;  % the must larger than cr  
    step = delta_c;
    c_right = beta_min;   %    right boundary of this function
    for m=2:N+1
        if c_right > beta(m)
            c_right = beta(m);
        end
    end
    c_temp = c;
    F0 = real(SECULAR_Ray(w,c_temp));
    c_temp  = c+step;
    F1 = real(SECULAR_Ray(w,c_temp)); 
   
    while (~(F0<0 && F1>0))
       
        if (c_temp > c_right)
            root_real = [];
            return
        else
            c_temp = c_temp +step;
            F0 = F1;
            F1 = real(SECULAR_Ray(w,c_temp));
         
        end
    end
   
    % c_temp must less than  c_love_max
    if (c_temp >= c_right) % must > = ,if not ,may cause x1 = x2，the answer is Inf
        c_temp = c_right;
        F1 = real(SECULAR_Ray(w,c_temp));
        if (F0/F1 > 0)  %Judge whether the interval has a solution
            root_real = [];
            return
        end    
        
    end
   
    % The interval of the root is found
    
    % find the root using secant method。
    m = 0;      
    n_iter = 500;   %Number of iterations
    x_min = c_temp-step;
    x1 = x_min;
    x2 = c_temp;
    while (m<n_iter)
        m = m+1;
        x_temp = x1-F0/(F1-F0)*(x2-x1);
        if (x1<x_min)
            x1 = x_min;
        elseif (x2<x_min)
            x2 = x_min;
        end
        if (abs(x_temp-x2) > Epsilon)
            x1 = x2;
            x2 = x_temp;
            F0 = real(SECULAR_Ray(w,x1)); 
            F1 = real(SECULAR_Ray(w,x2)); 
        else
            break
        end
    end

    if (abs(imag(SECULAR_Ray(w,x_temp)))<10^(-4))
       
        root_real =  x_temp;
    else
        root_real = [];
    end
    
end
    
    
    
    
    
    







function c_fund = c_fundamental_appr(w)
    %Calculation of approximate fundamental Rayleigh wave velocity
    global  beta alpha miu  N z_bottom  cr ;
    global nu gama ;
    c = cr ;  
    k=w/c;
    
    % secular function at single frequency and single velocity(R(c))
   
    %  Initialization parameters
    nu = zeros(1,N+1);       %  v,    (k^2-(w/beta)^2)^(1/2)
    gama = zeros(1,N+1);     % gama,  (k^2-(w/alpha)^2)^(1/2)
    X = zeros(1,N+1);        % Chi, X = k^2 + v^2 
    E = zeros(4,4,N+1);      % 3D array,3th dimension represents the number of layers
    E_t = zeros(4,4,N+1);    % E_t is a temporary matrix to get E matrix ; E=1/w*E_t（9）
    
    E11 = zeros(2,2,N+1);
    E12 = zeros(2,2,N+1);
    E21 = zeros(2,2,N+1);
    E22 = zeros(2,2,N+1);
    RT_matrix = zeros(4,4,N-1);  % R/T matrix 
    RT_N = zeros(4,2);           % j=N  (Td Rdu)T (24b)
    Td  = zeros(2,2,N);          % R/T  coefficients     (24a)
    Rud = zeros(2,2,N);
    Rdu = zeros(2,2,N);
    Tu  = zeros(2,2,N);
    Rdu_G = zeros(2,2,N);        % generalized reflection and transmission coefficients
    Td_G  = zeros(2,2,N); 
    Zero_M_2 = zeros(2,2);       % 2*2 ，0 matrix 
    I_2  =  diag([1,1]);           % Identity matrix
    
    
    
    % when c == beta(j) ,the secular function is singular .in this case, i plus
    % a small number to c 
    % i think it's useless in this function
   %{
    for n= 1:N+1
        if c == beta(n)
            c = c+Epsilon;
        end
    end
    %}
    
    for m = 1:N+1
        nu(m)   = (k^2-(w/beta(m))^2)^(1/2);
        gama(m) = (k^2-(w/alpha(m))^2)^(1/2);
        X(m) = k^2+nu(m)^2;   
        E_t(1,1,m) = alpha(m)*k;                 % (9)
        E_t(1,2,m) = beta(m)*nu(m);
        E_t(1,3,m) = E_t(1,1,m);
        E_t(1,4,m) = E_t(1,2,m);
        E_t(2,1,m) = alpha(m)*gama(m) ;
        E_t(2,2,m) = beta(m)*k;
        E_t(2,3,m) = -E_t(2,1,m);
        E_t(2,4,m) = -E_t(2,2,m);
        E_t(3,1,m) = -2*alpha(m)*miu(m)*k*gama(m);
        E_t(3,2,m) = -beta(m)*miu(m)*X(m);
        E_t(3,3,m) = -E_t(3,1,m) ;
        E_t(3,4,m) = -E_t(3,2,m) ;
        E_t(4,1,m) = -alpha(m)*miu(m)*X(m);
        E_t(4,2,m) = -2*beta(m)*miu(m)*k*nu(m);
        E_t(4,3,m) = E_t(4,1,m);
        E_t(4,4,m) = E_t(4,2,m);
        
    end
    E = 1/w*E_t ;   
    E11 = E(1:2,1:2,:);
    E12 = E(1:2,3:4,:);
    E21 = E(3:4,1:2,:);
    E22 = E(3:4,3:4,:);
    for m=1:N-1
        
        Lm_d_M =  Lambda_Ray_d(z_bottom(m),m);      % (24a)
        Lm_u_M =  Lambda_Ray_u(z_bottom(m),m+1);    % (24a)
        L_du_M = [Lm_d_M,Zero_M_2;Zero_M_2,Lm_u_M]; % (24a) item 3 
        E_inv_M = [E11(:,:,m+1),-E12(:,:,m);E21(:,:,m+1),-E22(:,:,m)]; % (24a) item 1
        E_M = [E11(:,:,m),-E12(:,:,m+1);E21(:,:,m),-E22(:,:,m+1)];     % (24a) item 2
        RT_matrix(:,:,m) =  E_inv_M\E_M*L_du_M ;           % (24a)
        
        Td(:,:,m)  =  RT_matrix(1:2,1:2,m);
       
        Rud(:,:,m) =  RT_matrix(1:2,3:4,m);
        Rdu(:,:,m) =  RT_matrix(3:4,1:2,m);
        Tu(:,:,m)  =  RT_matrix(3:4,3:4,m);
    end
    E_inv_M = [E11(:,:,N+1),-E12(:,:,N);E21(:,:,N+1),-E22(:,:,N)]; % (24b) item 1
    Lm_d_N  =  Lambda_Ray_d(z_bottom(N),N);   
    E_M =  [E11(:,:,N)*Lm_d_N;E21(:,:,N)*Lm_d_N];
    
    % j= N, interface
    RT_N = E_inv_M\E_M;                 %(24b)
    Td(:,:,N)  = RT_N(1:2,1:2);         %(24b)
    Rdu(:,:,N) = RT_N(3:4,1:2);         %(24b)
    
    % generalized R/T coefficients (Rdu_G Td_G)
    % j=N, interface
    Rdu_G(:,:,N) = Rdu(:,:,N);
    Td_G(:,:,N)  = Td(:,:,N);
    
    % j=1,2,3...N-1 interface
    for m = N-1:-1:1
        Td_G(:,:,m) = (I_2-Rud(:,:,m)*Rdu_G(:,:,m+1))\Td(:,:,m);  % （28） note! it is \ ,not / !!!
        Rdu_G(:,:,m)= Rdu(:,:,m)+Tu(:,:,m)*Rdu_G(:,:,m+1)*Td_G(:,:,m); % (28)
    end
    A  = E21(:,:,1)/w;    % A6
    lmu_0 = Lambda_Ray_u(0,1);
    B  = E22(:,:,1)*lmu_0*Rdu_G(:,:,1)/w;
    g_cr = A(1,1)*B(2,2)+A(2,2)*B(1,1)-A(1,2)*B(2,1)+det(B);
    R_deri_cr = R_deri(c);
    c_fund = c-(g_cr)/(miu(1)^2*alpha(1)*beta(1)*R_deri_cr);     % A9




end

function R_c = R(c)
    %    A1
    
    global beta alpha;
    a   = alpha(1);
    b   = beta(1);
    c_a = (1/c)^2-(1/a)^2;
    c_b = (1/c)^2-(1/b)^2;
    R_c = 4/c^2*c_a^(1/2)*c_b^(1/2)-(2/c^2-1/b^2)^2;
    

end
function R_deri_c = R_deri(c)
    %   A1   Derivative of R(c)
    global beta alpha;
    a   = alpha(1);
    b   = beta(1);  
    c_a = (1/c)^2-(1/a)^2;
    c_b = (1/c)^2-(1/b)^2;
    I_1 = -8/c^3*(c_a*c_b)^(1/2);   %  First term of derivative polynomial
    I_2 = -4/c^5*(c_a^(-1/2)*c_b^(1/2)+c_a^(1/2)*c_b^(-1/2));  %  second term
    I_3 = 8/c^3*(2/c^2-1/b^2);
    R_deri_c = I_1+I_2+I_3;


end




function  rayleigh_equation_value = Rayleigh_equation(c)
    % When Rayleigh_equation is equal to 0, the c is CR( Rayleigh velocity)
    global beta alpha;
    b   = beta(1);   %beta(1)
    a   = alpha(1);
    c_b = c/b;
    b_a = b/a;
    rayleigh_equation_value = c_b^6-8*c_b^4+c^2*(24/(b^2)-16/(a^2)) ...
        -16*(1-b_a^2);      % （5-2-18）introduction to seismology Wan Yongge
        
    

end

function CR = getCR()
    % get the Rayleigh velocity, cR, for a homogeneous half-space media with the 
    % same elastic properties as our top layer
    % Finding root by chord-truncating method
    global beta alpha  c_Ray_min c_Ray_max delta_c  Epsilon;
    c1  =  c_Ray_min;
    f1  =  Rayleigh_equation(c1);
    c2  =  c1+delta_c;
    f2  =  Rayleigh_equation(c2);
    while(f1/f2>0)
        if(c2>c_Ray_max)
            fprintf('the range of Rayleigh velocity is wrong');
            return
        end
        c1 = c2;
        c2 = c2 + delta_c;
        f1  =  Rayleigh_equation(c1);
        f2  =  Rayleigh_equation(c2);
    end
   m = 0;
   while( m<500)
       m = m+1;
       t = c1-f1*(c2-c1)/(f2-f1);
       if(abs(t-c2)< Epsilon)
           
           break
       else
           c1 = c2;
           c2 = t;
           f1  =  Rayleigh_equation(c1);
           f2  =  Rayleigh_equation(c2);
       end
   end
   CR = t;

end





 function  lambda_Ray_d = Lambda_Ray_d(z,j)     % j=1,2...N (22a)
    global nu gama z_bottom;
    lambda_Ray_d = zeros(2,2);
    if (j == 1)
        lmd = [exp(-gama(j)*z),exp(-nu(j)*z)];
        
    else
        lmd = [exp(-gama(j)*(z-z_bottom(j-1))),exp(-nu(j)*(z-z_bottom(j-1)))];
    end
    lambda_Ray_d = diag(lmd);

 end
 
 function  lambda_Ray_u = Lambda_Ray_u(z,j)     % j=1,2...N (22a)
    global nu gama z_bottom;
    lambda_Ray_u = zeros(2,2);
    
    lmu = [exp(-gama(j)*(z_bottom(j)-z)),exp(-nu(j)*(z_bottom(j)-z))];
   
    lambda_Ray_u = diag(lmu);

 end
 
 
 
 
 
 
function secular_Ray_w_c = SECULAR_Ray(w,c)  
    % secular function vlaue at (w,c) 
    global N z_bottom  beta miu  Epsilon alpha;
    % secular function at single frequency and single velocity
    global nu gama ;
    %  Initialization parameters
    nu = zeros(1,N+1);       %  v,    (k^2-(w/beta)^2)^(1/2)
    gama = zeros(1,N+1);     % gama,  (k^2-(w/alpha)^2)^(1/2)
    X = zeros(1,N+1);        % Chi, X = k^2 + v^2 
    E = zeros(4,4,N+1);      % 3D array,3th dimension represents the number of layers
    E_t = zeros(4,4,N+1);    % E_t is a temporary matrix to get E matrix ; E=1/w*E_t（9）
    
    E11 = zeros(2,2,N+1);
    E12 = zeros(2,2,N+1);
    E21 = zeros(2,2,N+1);
    E22 = zeros(2,2,N+1);
    RT_matrix = zeros(4,4,N-1);  % R/T matrix 
    RT_N = zeros(4,2);           % j=N  (Td Rdu)T (24b)
    Td  = zeros(2,2,N);          % R/T  coefficients     (24a)
    Rud = zeros(2,2,N);
    Rdu = zeros(2,2,N);
    Tu  = zeros(2,2,N);
    Rdu_G = zeros(2,2,N);        % generalized reflection and transmission coefficients
    Td_G  = zeros(2,2,N); 
    Zero_M_2 = zeros(2,2);       % 2*2 ，0 matrix 
    I_2  =  diag([1,1]);           % Identity matrix
    
   
    
    % when c == beta(j) ,the secular function is singular .in this case, i plus
    % a small number to c
    
    for n= 1:N+1
        if (abs(c-beta(n))<(Epsilon)||abs(c-alpha(n))<(Epsilon))
            
            c = c+Epsilon;
            
        end
    end

     k=w/c;
    for m = 1:N+1
        nu(m)   = (k^2-(w/beta(m))^2)^(1/2);
        gama(m) = (k^2-(w/alpha(m))^2)^(1/2);
        X(m) = k^2+nu(m)^2;   
        E_t(1,1,m) = alpha(m)*k;                 % (9)
        E_t(1,2,m) = beta(m)*nu(m);
        E_t(1,3,m) = E_t(1,1,m);
        E_t(1,4,m) = E_t(1,2,m);
        E_t(2,1,m) = alpha(m)*gama(m) ;
        E_t(2,2,m) = beta(m)*k;
        E_t(2,3,m) = -E_t(2,1,m);
        E_t(2,4,m) = -E_t(2,2,m);
        E_t(3,1,m) = -2*alpha(m)*miu(m)*k*gama(m);
        E_t(3,2,m) = -beta(m)*miu(m)*X(m);
        E_t(3,3,m) = -E_t(3,1,m) ;
        E_t(3,4,m) = -E_t(3,2,m) ;
        E_t(4,1,m) = -alpha(m)*miu(m)*X(m);
        E_t(4,2,m) = -2*beta(m)*miu(m)*k*nu(m);
        E_t(4,3,m) = E_t(4,1,m);
        E_t(4,4,m) = E_t(4,2,m);
        
    end
    E = 1/w*E_t;    
    E11 = E(1:2,1:2,:);
    E12 = E(1:2,3:4,:);
    E21 = E(3:4,1:2,:);
    E22 = E(3:4,3:4,:);
    for m=1:N-1
        
        Lm_d_M =  Lambda_Ray_d(z_bottom(m),m);      % (24a)
        Lm_u_M =  Lambda_Ray_u(z_bottom(m),m+1);    % (24a)
        L_du_M = [Lm_d_M,Zero_M_2;Zero_M_2,Lm_u_M]; % (24a) item 3 
        E_inv_M = [E11(:,:,m+1),-E12(:,:,m);E21(:,:,m+1),-E22(:,:,m)]; % (24a) item 1
        E_M = [E11(:,:,m),-E12(:,:,m+1);E21(:,:,m),-E22(:,:,m+1)];     % (24a) item 2
        RT_matrix(:,:,m) =  E_inv_M\E_M*L_du_M ;           % (24a)
        
        Td(:,:,m)  =  RT_matrix(1:2,1:2,m);
       
        Rud(:,:,m) =  RT_matrix(1:2,3:4,m);
        Rdu(:,:,m) =  RT_matrix(3:4,1:2,m);
        Tu(:,:,m)  =  RT_matrix(3:4,3:4,m);
    end
    E_inv_M = [E11(:,:,N+1),-E12(:,:,N);E21(:,:,N+1),-E22(:,:,N)];  % (24b) item 1
    Lm_d_N  =  Lambda_Ray_d(z_bottom(N),N);   
    E_M =  [E11(:,:,N)*Lm_d_N;E21(:,:,N)*Lm_d_N];

    % j= N, interface
    RT_N = E_inv_M\E_M ;                 %(24b)
    Td(:,:,N)  = RT_N(1:2,1:2);         %(24b)
    Rdu(:,:,N) = RT_N(3:4,1:2);         %(24b)
    
    % generalized R/T coefficients (Rdu_G Td_G)
    % j=N, interface
    Rdu_G(:,:,N) = Rdu(:,:,N);
    Td_G(:,:,N)  = Td(:,:,N);
    
    % j=1,2,3...N-1 interface
    for m = N-1:-1:1
        Td_G(:,:,m) = (I_2-Rud(:,:,m)*Rdu_G(:,:,m+1))\Td(:,:,m);  % （28） note! it is \ ,not / !!!
        Rdu_G(:,:,m)= Rdu(:,:,m)+Tu(:,:,m)*Rdu_G(:,:,m+1)*Td_G(:,:,m); % (28)
    end
    % j=0   Rud_G(0)   (27)
    Lm_u_0 =  Lambda_Ray_u(0,1);    % (24a)
    Rud_G0  = -E21(:,:,1)\E22(:,:,1)*Lm_u_0;
    
   
    % secular function
    secular_Ray_w_c = det(I_2-Rud_G0*Rdu_G(:,:,1));      %      (35)
    
    
    
    
    
    
end

function [root_one,x_right]=findroot_1(w0,c_start,step) 
% x_left is the right boundary of root interval
% case : w = w0,find a sinale root
    global  c_Ray_max Epsilon;
    
    %  Find the interval of root
    F0 = imag(SECULAR_Ray(w0,c_start));
    c_temp = c_start + step;
    F1 = imag(SECULAR_Ray(w0,c_temp)); 
    
    while (F0/F1 > 0)
       
        if (c_temp > c_Ray_max)
            root_one = [];
            x_right = [];
            return
        else
            c_temp = c_temp +step;
            F0 = F1;
            F1 = imag(SECULAR_Ray(w0,c_temp));

        end
    end
    x_right = c_temp;
    % c_temp must less than  c_love_max
    if (c_temp >= c_Ray_max) % must > = ,if not ,may cause x1 = x2，the answer is Inf
        c_temp = c_Ray_max;
        F1 = imag(SECULAR_Ray(w0,c_temp));
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
            F0 = imag(SECULAR_Ray(w0,x1)); 
            F1 = imag(SECULAR_Ray(w0,x2)); 
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
        F = imag(SECULAR_Ray(w0,c_temp)); 
        if abs(F) > 0.01
            break
        end
        c_temp =  c_temp+step;
    end
    c_begin = c_temp; 
end

function root_all = findroot_all(w0)
    % find all roots when w = w0
    global c_Ray_max  delta_c delta_samll_c Epsilon  ;
    root_all = [];
    %c_start = find_begin(w0);
    c_start = find_begin(w0);
    % first root
    step = delta_samll_c;
    
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
         if (c_right>= c_Ray_max)
             if real(SECULAR_Ray(w0,c_temp)) < 0.3    % judge real part
                 root_all = [root_all,c_temp];  
             end
             
             return
         end
         c_start = c_right;
         if real(SECULAR_Ray(w0,c_temp)) < 0.3    % judge real part
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
         if (c_right>= c_Ray_max)
             if real(SECULAR_Ray(w0,c_temp)) < 0.3    % judge real part
                 root_all = [root_all,c_temp];  
             end
             
             return
         end
         c_start = c_right;
         if real(SECULAR_Ray(w0,c_temp)) < 0.3    % judge real part
              root_all = [root_all,c_temp];
  
         end
    end
end




