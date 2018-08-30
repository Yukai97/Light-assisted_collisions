clear;
clc;

t = -10*pi:pi/250:10*pi;
x1 = (cos(2*t).^2).*sin(t);
y1 = (sin(2*t).^2).*cos(t);
x2 = (sin(2*t).^2).*cos(t);
y2 = (cos(2*t).^2).*sin(t);
L = length(t);

for i = 1:L
    plot3(x1(i),y1(i),t(i),'bo','MarkerSize',10,'MarkerFaceColor','blue');
    hold on;
    plot3(x2(i),y2(i),t(i),'ro','MarkerSize',10,'MarkerFaceColor','red');
    hold off;
    axis([-1 1 -1 1 -10*pi 10*pi]);
    grid on;
    drawnow limitrate
end