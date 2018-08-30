function [] = plot_trajectories(rvec1,rvec2,N_t,spacing,w0,z_R)

rvecp1 = [rvec1(1:spacing:N_t,:);rvec1(N_t,:)];
rvecp2 = [rvec2(1:spacing:N_t,:);rvec1(N_t,:)];
L = length(rvecp1);

for i = 1:L
    plot3(rvecp1(i,1),rvecp1(i,2),rvecp1(i,3),'bo','MarkerSize',10,'MarkerFaceColor','blue');
    hold on;
    plot3(rvecp2(i,1),rvecp2(i,2),rvecp2(i,3),'ro','MarkerSize',10,'MarkerFaceColor','red');
    hold off;
    grid on;
    %axis([-w0,w0,-w0,w0,-z_R,z_R]);
    drawnow;
end




end

