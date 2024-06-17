function coordinates = localCoordinateSystem(prox,dist,planar)

x_static = prox - dist;

temp1 = planar - dist;

y_static = cross(temp1,x_static);

z_static = cross(x_static, y_static);

distx = norm(x_static);
disty = norm(y_static);
distz = norm(z_static);

x_static = x_static/distx;
y_static = y_static/disty;
z_static = z_static/distz;

coordinates = horzcat(x_static', y_static', z_static');

end