function [X,Y,N] =Runge_Kutta(f,h,a,b,y0)
%利用四阶Runge-Kutta公式求解常微分方程
%f为函数，h为步长,a为区间左端点，b为区间右端点,y0为初值
N=floor((b-a)/h+1.001);     %此处加1.001再取整，是防止(b-a)/h出现X.9999情况，实际应取X+1，但直接用floor将取X，故此处加上1.001总体取整
X=zeros(1,N);
Y=zeros(1,N);

for i=1:1:N
    X(i)=a+(i-1)*h;
end                             %确定分点

Y(1)=y0;
for i=1:1:N-1
    k1=f(X(i),Y(i));
    k2=f(X(i)+h/2,Y(i)+h/2*k1);
    k3=f(X(i)+h/2,Y(i)+h/2*k2);
    k4=f(X(i)+h,Y(i)+h*k3);
    Y(i+1)=Y(i)+h/6*(k1+2*k2+2*k3+k4);       %由四阶Runge-Kutta公式求解微分方程
end

end

