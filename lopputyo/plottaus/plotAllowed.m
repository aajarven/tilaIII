function plotAllowed(mu, C)
    v = -2:0.005:2;
    [x,y] = meshgrid(v);
    ineq = 0.5*(x.^2+y.^2) + (1-mu)./(sqrt((mu+x).^2+y.^2)) + mu./sqrt((x-1+mu).^2+y.^2) >= C/2;
    f = double(ineq);
    s = surf(x,y,f);
    hold on;
    view(0,90);
    s.EdgeColor = 'none';
    pbaspect([1 1 1]);
    colormap gray
    title({'Forbidden areas (in black) for orbits', strcat(strcat('having µ=', num2str(mu)), strcat(' and C=', num2str(C)))});
    xlabel('x')
    ylabel('y')
    caxis([0 1])
    scatter3(-mu, 0, 1, 50, 'o', 'filled', 'y', 'MarkerEdgeColor', 'k');
    scatter3(1-mu, 0, 1, 20, 'o', 'filled', 'r');
    hold off;
end