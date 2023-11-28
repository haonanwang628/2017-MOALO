function [ obj ] =test_function(POP,c,m,F)  %基准测试问题
% Usage: [ obj ] =compute_objectives(POP,c,problem_name)
%
% Input:
% problem_name  - Benchmark MOP (ZDT, DTLZ, WFG Problems)
% c             -No. of Decision Variables     %决策变量序列号
% POP           -Population of Decision Variables   %决策变量种群
% m             %目标维数
% Output: 
% obj           - Calculated Objective Values  %目标函数值
%
    %%%%    Authors:    Handing Wang, Licheng Jiao, Xin Yao
    %%%%    Xidian University, China, and University of Birmingham, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com, X.Yao@cs.bham.ac.uk
    %%%%    WEBSITE:    http://www.cs.bham.ac.uk/~xin/
    %%%%    DATE:       August 2014
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Licheng Jiao, Xin Yao, An Improved Two-Archive Algorithm for Many-Objective Optimization, Evolutionary Computation, IEEE Transactions on, Accepted, 10.1109/TEVC.2014.2350987.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
obj=[];
n=size(POP,1);
switch   F         %problem_name
    case 'F1'       %ZDT1
        obj(:,1)=POP(:,1);
        g=ones(n,1)+9*sum(POP(:,2:c),2)/(c-1);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^0.5);
   
    case 'F2'       %ZDT2
        obj(:,1)=POP(:,1);
        g=ones(n,1)+9*sum(POP(:,2:c),2)/(c-1);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^2);
   
    case 'F3'       %ZDT3
        obj(:,1)=POP(:,1);
        g=ones(n,1)+9*sum(POP(:,2:c),2)/(c-1);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^0.5-(POP(:,1)./g).*(sin(10*pi*POP(:,1))));
    
    case 'F4'     %ZDT4
        obj(:,1)=POP(:,1);
        g=ones(n,1)+10*(c-1)*ones(n,1)+sum((POP(:,2:c).^2-10*cos(4*pi*POP(:,2:c))),2);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^0.5);
        
    case 'F5'     %ZDT5
        u = zeros(size(POP,1),1+(size(POP,2)-30)/5);
        u(:,1) = sum(POP(:,1:30),2);
        for i = 2 : size(u,2)
            u(:,i) = sum(POP(:,(i-2)*5+31:(i-2)*5+35),2);
        end
        v           = zeros(size(u));
        v(u<5)      = 2 + u(u<5);
        v(u==5)     = 1;
        obj(:,1) = 1 + u(:,1);
        g           = sum(v(:,2:c),2);
        h           = 1./obj(:,1);
        obj(:,2) = g.*h;
        
    case 'F6'      %ZDT6
        obj(:,1)=ones(n,1)-exp(-4*POP(:,1)).*((sin(6*pi*POP(:,1))).^6);
        g=ones(n,1)+9*((sum(POP(:,2:c),2)/(c-1)).^0.25);
        obj(:,2)=g.*(ones(n,1)-(obj(:,1)./g).^2);

    case 'F7'   %DTLZ1_m          
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=0.5*prod(POP(:,1:m-1),2).*(1+g);        %POP中前m-1列每行连乘
        for i=2:m-1
            obj(:,i)=0.5*prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=0.5*(1-POP(:,1)).*(1+g);

    case 'F8'     %DTLZ2_m
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);

    case 'F9'     %DTLZ3_m
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
   
     case 'F10'   %DTLZ4_m
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*(POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*(POP(:,1:m-i).^100)),2).*(sin(0.5*pi*(POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*(POP(:,1).^100))).*(1+g);
       
    case 'F11'   %DTLZ5_m
%         g=sum((POP(:,m:c)-0.5).^2,2);
%         Temp=repmat(g,1,m-2);
%         theta(:,2:m-1) = pi*(1+2*Temp.*POP(:,2:m-1))./(4*(1+Temp));
%         obj(:,1)=prod(cos(0.5*pi*theta(:,1:m-1)),2).*(1+g);
%         for i=2:m-1
%             obj(:,i)=prod(cos(0.5*pi*theta(:,1:m-i)),2).*(sin(0.5*pi*theta(:,m+1-i))).*(1+g);
%         end
%         obj(:,m)=(sin(0.5*pi*theta(:,1))).*(1+g); 
         g      = sum((POP(:,m:end)-0.5).^2,2);
         Temp   = repmat(g,1,m-2);
         POP(:,2:m-1) = pi*(1+2*Temp.*POP(:,2:m-1))./(4*(1+Temp));
         obj = repmat(1+g,1,m).*fliplr(cumprod([ones(size(g,1),1),cos(POP(:,1:m-1)*pi/2)],2)).*[ones(size(g,1),1),sin(POP(:,m-1:-1:1)*pi/2)];
        
    case 'F12'   %DTLZ6_m
         g = sum(POP(:,m:c).^0.1,2);
         Temp = repmat(g,1,m-2);
         POP(:,2:m-1) = pi*(1+2*Temp.*POP(:,2:m-1))./(4*(1+Temp));
         obj = repmat(1+g,1,m).*fliplr(cumprod([ones(size(g,1),1),cos(POP(:,1:m-1)*pi/2)],2)).*[ones(size(g,1),1),sin(POP(:,m-1:-1:1)*pi/2)];
        
        
    case 'F13'    %DTLZ7_m
         obj = zeros(size(POP,1),m);
         g = 1+9*mean(POP(:,m:c),2);
         obj(:,1:m-1) = POP(:,1:m-1);
         obj(:,m) = (1+g).*(m-sum(obj(:,1:m-1)./(1+repmat(g,1,m-1)).*(1+sin(3*pi.*obj(:,1:m-1))),2));
    
    case 'F14'    %WFG1_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c1 = 1;
            S = 2 : 2 : 2*m;
            A = ones(1,m-1);

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);

            t1 = zeros(N,K+L);
            t1(:,1:K)     = z01(:,1:K);
            t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

            t2 = zeros(N,K+L);
            t2(:,1:K)     = t1(:,1:K);
            t2(:,K+1:end) = b_flat(t1(:,K+1:end),0.8,0.75,0.85);

            t3 = zeros(N,K+L);
            t3 = b_poly(t2,0.02);

            t4 = zeros(N,m);
            for i = 1 : m-1
                t4(:,i) = r_sum(t3(:,(i-1)*K/(m-1)+1:i*K/(m-1)),2*((i-1)*K/(m-1)+1):2:2*i*K/(m-1));
            end
            t4(:,m) = r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t4(:,m),A(i)).*(t4(:,i)-0.5)+0.5;
            end
            x(:,m) = t4(:,m);

            h      = convex(x);
            h(:,m) = mixed(x);
            obj = repmat(c1*x(:,m),1,m) + repmat(S,N,1).*h;
    
    case 'F15'    %WFG2_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c1 = 1;
            S = 2 : 2 : 2*m;
            A = ones(1,m-1);

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1(:,1:K)     = z01(:,1:K);
            t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

            t2 = zeros(N,K+L/2);
            t2(:,1:K) = t1(:,1:K);
            % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
            t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
            % ---------------------------------------------------------
            
            t3 = zeros(N,m);
            for i = 1 : m-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(m-1)+1:i*K/(m-1)),ones(1,K/(m-1)));
            end
            t3(:,m) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t3(:,m),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,m) = t3(:,m);

            h      = convex(x);
            h(:,m) = disc(x);
            obj = repmat(c1*x(:,m),1,m) + repmat(S,N,1).*h;
        
    case 'F16'     %WFG3_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c1 = 1;
            S = 2 : 2 : 2*m;
            A = [1,zeros(1,m-2)];

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);

            t1 = zeros(N,K+L);
            t1(:,1:K)     = z01(:,1:K);
            t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

            t2 = zeros(N,K+L/2);
            t2(:,1:K) = t1(:,1:K);
            % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
            t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
            % ---------------------------------------------------------
            
            t3 = zeros(N,m);
            for i = 1 : m-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(m-1)+1:i*K/(m-1)),ones(1,K/(m-1)));
            end
            t3(:,m) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t3(:,m),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,m) = t3(:,m);

            h      = linear(x);
            obj = repmat(c1*x(:,m),1,m) + repmat(S,N,1).*h;
        
    case 'F17'     %WFG4_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c = 1;
            S = 2 : 2 : 2*m;
            A = ones(1,m-1);

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);

            t1 = zeros(N,K+L);
            t1 = s_multi(z01,30,10,0.35);

            t2 = zeros(N,m);
            for i = 1 : m-1
                t2(:,i) = r_sum(t1(:,(i-1)*K/(m-1)+1:i*K/(m-1)),ones(1,K/(m-1)));
            end
            t2(:,m) = r_sum(t1(:,K+1:K+L),ones(1,L));

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t2(:,m),A(:,i)).*(t2(:,i)-0.5)+0.5;
            end
            x(:,m) = t2(:,m);

            h = concave(x);
            obj = repmat(c*x(:,m),1,m) + repmat(S,N,1).*h;
        
    case 'F18'      %WFG5_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c1 = 1;
            S = 2 : 2 : 2*m;
            A = ones(1,m-1);

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1 = s_decept(z01,0.35,0.001,0.05);

            t2 = zeros(N,m);
            for i = 1 : m-1
                t2(:,i) = r_sum(t1(:,(i-1)*K/(m-1)+1:i*K/(m-1)),ones(1,K/(m-1)));
            end
            t2(:,m) = r_sum(t1(:,K+1:K+L),ones(1,L));

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t2(:,m),A(:,i)).*(t2(:,i)-0.5)+0.5;
            end
            x(:,m) = t2(:,m);

            h = concave(x);
            obj = repmat(c1*x(:,m),1,m) + repmat(S,N,1).*h;
        
    case 'F19'      %WFG6_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c1 = 1;
            S = 2 : 2 : 2*m;
            A = ones(1,m-1);

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1(:,1:K)     = z01(:,1:K);
            t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

            t2 = zeros(N,m);
            for i = 1 : m-1
                t2(:,i) = r_nonsep(t1(:,(i-1)*K/(m-1)+1:i*K/(m-1)),K/(m-1));
            end
            % Same as <t2(:,M)=r_nonsep(t1(:,K+1:end),L)>
            SUM = zeros(N,1);
            for i = K+1 : K+L-1
                for j = i+1 : K+L
                    SUM = SUM + abs(t1(:,i)-t1(:,j));
                end
            end
            t2(:,m) = (sum(t1(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
            % -------------------------------------------

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t2(:,m),A(:,i)).*(t2(:,i)-0.5)+0.5;
            end
            x(:,m) = t2(:,m);

            h = concave(x);
            obj = repmat(c1*x(:,m),1,m) + repmat(S,N,1).*h;
        
    case 'F20'       %WFG7_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c1 = 1;
            S = 2 : 2 : 2*m;
            A = ones(1,m-1);

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
            Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
            t1(:,1:K) = z01(:,1:K).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K)).*abs(floor(0.5-Y(:,1:K))+0.98/49.98)));
            % ------------------------------------------------------------------------------------------
            t1(:,K+1:end) = z01(:,K+1:end);

            t2 = zeros(N,K+L);
            t2(:,1:K)     = t1(:,1:K);
            t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);

            t3 = zeros(N,m);
            for i = 1 : m-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(m-1)+1:i*K/(m-1)),ones(1,K/(m-1)));
            end
            t3(:,m) = r_sum(t2(:,K+1:K+L),ones(1,L));

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t3(:,m),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,m) = t3(:,m);

            h = concave(x);
            obj = repmat(c1*x(:,m),1,m) + repmat(S,N,1).*h;
        
    case 'F21'       %WFG8_m
            N = size(POP,1);
            K = m-1;
            L = c - K;
            c1 = 1;
            S = 2 : 2 : 2*m;
            A = ones(1,m-1);

            z01 = POP./repmat(2:2:size(POP,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1(:,1:K) = z01(:,1:K);
            % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,1:i-1),ones(1,i-1)),0.98/49.98,0.02,50)>
            Y = (cumsum(z01,2)-z01)./repmat(0:K+L-1,N,1);
            t1(:,K+1:K+L) = z01(:,K+1:K+L).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,K+1:K+L)).*abs(floor(0.5-Y(:,K+1:K+L))+0.98/49.98))); 
            % --------------------------------------------------------------------------------------

            t2 = zeros(N,K+L);
            t2(:,1:K)     = t1(:,1:K);
            t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);

            t3 = zeros(N,m);
            for i = 1 : m-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(m-1)+1:i*K/(m-1)),ones(1,K/(m-1)));
            end
            t3(:,m) = r_sum(t2(:,K+1:K+L),ones(1,L));

            x = zeros(N,m);
            for i = 1 : m-1
                x(:,i) = max(t3(:,m),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,m) = t3(:,m);

            h = concave(x);
            obj = repmat(c1*x(:,m),1,m) + repmat(S,N,1).*h;
    
    
    
end


end

function Output = b_param(y,Y,A,B,C)
    Output = y.^(B+(C-B)*(A-(1-2*Y).*abs(floor(0.5-Y)+A)));
end

function Output = s_linear(y,A)
    Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = b_flat(y,A,B,C)
    Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
    Output = roundn(Output,-6);
end

function Output = b_poly(y,a)
    Output = y.^a;
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
    Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end

function Output = concave(x)
    Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end

function Output = r_nonsep(y,A)
    Output = zeros(size(y,1),1);
    for j = 1 : size(y,2)
        Temp = zeros(size(y,1),1);
        for k = 0 : A-2
            Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
        end
        Output = Output+y(:,j)+Temp;
    end
    Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end

function Output = s_decept(y,A,B,C)
    Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end

function Output = s_multi(y,A,B,C)
    Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end

function Output = disc(x)
    Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end