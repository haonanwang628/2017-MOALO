function  true_obj = TPF(m,N,F)            %’Ê µPareto«∞—ÿ
 switch F
     case 'F1'            %ZDT1
          true_obj(:,1) = (0:1/(N-1):1)';
          true_obj(:,2) = 1 - true_obj(:,1).^0.5;
           
     case 'F2'           %ZDT2
          true_obj(:,1) = (0:1/(N-1):1)';
          true_obj(:,2) = 1 - true_obj(:,1).^2;
          
     case 'F3'           %ZDT3
          true_obj(:,1) = (0:1/(N-1):1)';
          true_obj(:,2) = 1 - true_obj(:,1).^0.5 - true_obj(:,1).*sin(10*pi*true_obj(:,1));
          true_obj      = true_obj(NDSort(true_obj,1)==1,:);
          
     case 'F4'           %ZDT4
          true_obj(:,1) = (0:1/(N-1):1)';
          true_obj(:,2) = 1 - true_obj(:,1).^0.5;
          
          
     case 'F5'        %ZDT5
          c = 80;
          true_obj(:,1) = 1 : 31;
          true_obj(:,2) = (c-30)./5./true_obj(:,1); 
         
     case  'F6'          %ZDT6
           minf1 = 0.280775;
           true_obj(:,1) = (minf1:(1-minf1)/(N-1):1)';
           true_obj(:,2) = 1 - true_obj(:,1).^2;
    
     case 'F7'        %DTLZ1_m
%          a = repmat((0:1/(N-1):1)',1,m-1);
%          b = ones(N,c-m+1)*0.5;
%          x = [a,b];
%          g=100*((c-m+1)*ones(N,1)+sum((x(:,m:c)-0.5).^2-cos(20*pi*(x(:,m:c)-0.5)),2));
%          true_obj(:,1)=0.5*prod(x(:,1:m-1),2).*(1+g);
%          for i=2:m-1
%             true_obj(:,i)=0.5*prod(x(:,1:m-i),2).*(1-x(:,m+1-i)).*(1+g);
%          end
%         true_obj(:,m)=0.5*(1-x(:,1)).*(1+g);
          true_obj = UniformPoint(N,m)/2;
        
      case 'F8'     %DTLZ2_m
%          a = repmat((0:1/(N-1):1)',1,m-1);
%          b = ones(N,c-m+1)*0.5;
%          x = [a,b];
%          g=sum((x(:,m:c)-0.5).^2,2);
%          true_obj(:,1)=prod(cos(0.5*pi*x(:,1:m-1)),2).*(1+g);
%          for i=2:m-1
%              true_obj(:,i)=prod(cos(0.5*pi*x(:,1:m-i)),2).*(sin(0.5*pi*x(:,m+1-i))).*(1+g);
%          end
%          true_obj(:,m)=(sin(0.5*pi*x(:,1))).*(1+g);
           true_obj = UniformPoint(N,m); 
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
         
     case 'F9'      %DTLZ3_m
%          a = repmat((0:1/(N-1):1)',1,m-1);
%          b = ones(N,c-m+1)*0.5;
%          x = [a,b];
%          g=100*((c-m+1)*ones(N,1)+sum((x(:,m:c)-0.5).^2-cos(20*pi*(x(:,m:c)-0.5)),2));
%          true_obj(:,1)=prod(cos(0.5*pi*x(:,1:m-1)),2).*(1+g);
%          for i=2:m-1
%              true_obj(:,i)=prod(cos(0.5*pi*x(:,1:m-i)),2).*(sin(0.5*pi*x(:,m+1-i))).*(1+g);
%          end
%          true_obj(:,m)=(sin(0.5*pi*x(:,1))).*(1+g);
           true_obj = UniformPoint(N,m);
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
        
      case 'F10'   %DTLZ4_m
%          a = repmat((0:1/(N-1):1)',1,m-1);
%          b = ones(N,c-m+1)*0.5;
%          x = [a,b];
%          g = sum((x(:,m:c)-0.5).^2,2);
%          true_obj(:,1)=prod(cos(0.5*pi*(x(:,1:m-1).^100)),2).*(1+g);
%          for i=2:m-1
%              true_obj(:,i)=prod(cos(0.5*pi*(x(:,1:m-i).^100)),2).*(sin(0.5*pi*(x(:,m+1-i).^100))).*(1+g);
%          end
%          true_obj(:,m)=(sin(0.5*pi*(x(:,1).^100))).*(1+g);
%            a = repmat((0:1/(N-1):1)',1,m-1);
%            true_obj(:,1:m-1)=a;
%            true_obj(:,m)=sqrt(1-sum((true_obj(:,1:m-1).^2),2));
           true_obj = UniformPoint(N,m);
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
           
           
       case 'F11'   %DTLZ5_m
%            a = repmat((0:1/(N-1):1)',1,m-1);
%            true_obj(:,1:m-1)=a;
%            true_obj(:,m)=sqrt(1-sum((true_obj(:,1:m-1).^2),2));
           true_obj = [0:1/(N-1):1;1:-1/(N-1):0]';
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,size(true_obj,2));
           true_obj = [true_obj(:,ones(1,m-2)),true_obj];
           true_obj = true_obj./sqrt(2).^repmat([m-2,m-2:-1:0],size(true_obj,1),1);
          
     case 'F12'   %DTLZ6_m
           true_obj = [0:1/(N-1):1;1:-1/(N-1):0]';
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,size(true_obj,2));
           true_obj = [true_obj(:,ones(1,m-2)),true_obj];
           true_obj = true_obj./sqrt(2).^repmat([m-2,m-2:-1:0],size(true_obj,1),1);
           
     case 'F13'     %DTLZ7_m
           interval     = [0,0.251412,0.631627,0.859401];
           median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
           X            = ReplicatePoint(N,m-1);
           X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
           X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
           true_obj     = [X,2*(m-sum(X/2.*(1+sin(3*pi.*X)),2))];
     
     case 'F14'     %WFG1_m
            true_obj = UniformPoint(N,m);
            c = ones(size(true_obj,1),m);
            for i = 1 : size(true_obj,1) 
                for j = 2 : m
                    temp = true_obj(i,j)/true_obj(i,1)*prod(1-c(i,m-j+2:m-1));
                    c(i,m-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                end
            end
            x = acos(c)*2/pi;
            temp = (1-sin(pi/2*x(:,2))).*true_obj(:,m)./true_obj(:,m-1);
            a = 0 : 0.0001 : 1;
            E = abs(temp*(1-cos(pi/2*a))-1+repmat(a+cos(10*pi*a+pi/2)/10/pi,size(x,1),1));
            [~,rank] = sort(E,2);
            for i = 1 : size(x,1)
                x(i,1) = a(min(rank(i,1:10)));
            end
            true_obj      = convex(x);
            true_obj(:,m) = mixed(x);
            true_obj      = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
           
     case 'F15'   %WFG2_m
            true_obj = UniformPoint(N,m);
            c = ones(size(true_obj,1),m);
            for i = 1 : size(true_obj,1) 
                for j = 2 : m
                    temp = true_obj(i,j)/true_obj(i,1)*prod(1-c(i,m-j+2:m-1));
                    c(i,m-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                end
            end
            x = acos(c)*2/pi;
            temp = (1-sin(pi/2*x(:,2))).*true_obj(:,m)./true_obj(:,m-1);
            a = 0 : 0.0001 : 1;
            E = abs(temp*(1-cos(pi/2*a))-1+repmat(a.*cos(5*pi*a).^2,size(x,1),1));
            [~,rank] = sort(E,2);
            for i = 1 : size(x,1)
                x(i,1) = a(min(rank(i,1:10)));
            end
            true_obj      = convex(x);
            true_obj(:,m) = disc(x);
            true_obj      = true_obj(NDSort(true_obj,1)==1,:);
            true_obj      = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
            
     case 'F16'   %WFG3_m
           X = (0:1/(N-1):1)';
           X = [X,zeros(N,m-2)+0.5,zeros(N,1)];
           true_obj = linear(X);
           true_obj = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
           
     case 'F17'   %WFG4_m
           true_obj = UniformPoint(N,m);
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
           true_obj = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
           
     case 'F18'   %WFG5_m
           true_obj = UniformPoint(N,m);
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
           true_obj = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
           
     case 'F19'   %WFG6_m
           true_obj = UniformPoint(N,m);
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
           true_obj = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
           
     case 'F20'   %WFG7_m
           true_obj = UniformPoint(N,m);
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
           true_obj = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
           
     case 'F21'    %WFG8_m
           true_obj = UniformPoint(N,m);
           true_obj = true_obj./repmat(sqrt(sum(true_obj.^2,2)),1,m);
           true_obj = repmat(2:2:2*m,size(true_obj,1),1).*true_obj;
         
 end
end


function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end


function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
    Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end

function Output = disc(x)
    Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end