function plotPotential(mu)
    v = -1.6:0.1:1.6;
    [x,y] = meshgrid(v);
    U = -0.5*(x.^2+y.^2) - (1-mu)./(sqrt((mu+x).^2+y.^2)) - mu./sqrt((x-1+mu).^2+y.^2);
    s = surf(x,y,U);
    hold on;
    view(20,30);
    %s.EdgeColor = 'none';
    pbaspect([1 1 1]);
    colormap summer
    xlabel('x')
    ylabel('y')
    zlabel('U')
    zlim([-5 0])
    caxis([-5, 0])
    hold off;
end