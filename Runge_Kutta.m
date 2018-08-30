function [X,Y,N] =Runge_Kutta(f,h,a,b,y0)
%�����Ľ�Runge-Kutta��ʽ��ⳣ΢�ַ���
%fΪ������hΪ����,aΪ������˵㣬bΪ�����Ҷ˵�,y0Ϊ��ֵ
N=floor((b-a)/h+1.001);     %�˴���1.001��ȡ�����Ƿ�ֹ(b-a)/h����X.9999�����ʵ��ӦȡX+1����ֱ����floor��ȡX���ʴ˴�����1.001����ȡ��
X=zeros(1,N);
Y=zeros(1,N);

for i=1:1:N
    X(i)=a+(i-1)*h;
end                             %ȷ���ֵ�

Y(1)=y0;
for i=1:1:N-1
    k1=f(X(i),Y(i));
    k2=f(X(i)+h/2,Y(i)+h/2*k1);
    k3=f(X(i)+h/2,Y(i)+h/2*k2);
    k4=f(X(i)+h,Y(i)+h*k3);
    Y(i+1)=Y(i)+h/6*(k1+2*k2+2*k3+k4);       %���Ľ�Runge-Kutta��ʽ���΢�ַ���
end

end

