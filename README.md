%least square method
clc
clear
x=[-2 0 2 4];
y=[17.5 42.7 102 139.5];
n=4;
sum1=0;
sum2=0;
sum3=0;
sum4=0;
for i=1:n
    sum1=sum1+x(i);
    sum2=sum2+y(i);
    sum3=sum3+(x(i)*y(i));
    sum4=sum4+(x(i)*x(i));
end
syms a;
syms b;
A=[n sum1; sum1 sum4];
X=[a;b];
B=[sum2 ; sum3];
AX=B;
sum2=(n*a)+(b*sum1);
sum3=(a*sum1)+(b*sum4);

X= inv(A)*B




xxxxxxxxxxxxxx



power method 

clc 
clear all 
x0=[1;-1;2] 
A=[2 1 1;1 2 1;1 1 2] 
tol=10^-3; 
N=1000; 
y=A*x0; 
m1=max(y); 
i=0; 
while i<=N 
i=i+1;  
x=(1/m1).*y; 
y=A*x; 
m2=max(y); 
if abs(m2-m1)<tol 
fprintf(' max eigenvalue =%f',m2); 
break; 
end 
m1=m2; 
end 
x



xxxxxxxxxx

clc
clear
n=input('number of data points');
for i=1:n
    x(i)=input('enter the value of x(i)');
    f(i)=input('enter the value of function f(i)');
end
b=input('enter the point of interpolation b');
for j=1:n
    d(1,j)=f(j);
end
for i=2:n
    for j=1:n-i+1
        d(i,j)=(d(i-1,j+1)-d(i-1,j))/(x(i-1+j)-x(j));
    end
end
sum=d(1,1);
product=1;
for i=2:n
 product=product*(b-x(i-1));
        sum=sum+d(i,1)*product;
end

fprintf('value of f(%f) is =%f',b,sum)





#gausselinination

clc;

clear all; 

a= [1,1,1 1;4,3,-1 6;3 ,5,3 4] 

n=size(a,1);

l=eye(n);

for j=1: n-1 

for i=j+1: n 

l(i,j)=a(i,j) / a(j,j); 

a(i,:)=a(i, :) - l(i,j) *a(j,:);

end

end

x(n)=a(n , n+1)/a(n,n);

for i = n-1: -1 :1 

sum=0; 

for j=i+1 : n 

sum=sum + a(i,j) *x(j); 

end 

x(i)=(a(i, n+1)-sum)/a(i , i) ;

end

x




gauss seidal


% Gauss-Seidel Method
% Solves Ax = b
clc; clear;

A = [4 -1 0;
    -1 4 -1;
     0 -1 3];      % Coefficient matrix
b = [15; 10; 10];  % RHS vector
x0 = [0; 0; 0];    % Initial guess

tol = 1e-6;        % Tolerance
maxIter = 100;     % Maximum iterations

n = length(b);
x = x0;

for k = 1:maxIter
    x_old = x;
    
    for i = 1:n
        % sum1 = a(i,1:i-1)*x(1:i-1)
        sum1 = 0;
        for j = 1:i-1
            sum1 = sum1 + A(i,j)*x(j);
        end
        
        % sum2 = a(i,i+1:n)*x_old(i+1:n)
        sum2 = 0;
        for j = i+1:n
            sum2 = sum2 + A(i,j)*x_old(j);
        end
        
        x(i) = (b(i) - sum1 - sum2)/A(i,i);
    end
    
    % Check convergence
    if norm(x - x_old, inf) < tol
        fprintf('Converged in %d iterations\n', k);
        break;
    end
end

disp('Solution vector x = ');
disp(x);

