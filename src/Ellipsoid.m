% test the delta and hevaside function by calculating
% ellipsoid volume and surface area

% some typical value of a,b,c and the corresponding reduced volume
%	a		b		c 		t_rv	n_rv
%	210		210		114		0.900	0.902
%	200		150		71		0.800	0.794
%	215		215		70		0.700	0.704
%	215		215		62.3	0.650	0.654
%	215		215		55.5	0.600	0.605
%	215		215		49.3	0.550	0.556
%	215		215		43.6	0.500	0.508

% ellipsoid axis a>=b>=c
a = 215;
b = 215;
c = 215;

Vol = 4 * pi * a * b * c / 3; % volume of ellipsoid

phi = atan2(sqrt(a^2-c^2),c);
cos_phi = c/a;
sin_phi = sqrt(a^2-c^2)/a;
k = a * sqrt(b^2-c^2) / b / sqrt(a^2-c^2);

Area = 2*pi*c^2 + 2*pi*a*b*(ellipticE(phi,k)*sin_phi^2 + ellipticF(phi,k)*cos_phi^2) / sin_phi;

Red_Vol = (3*Vol/4/pi) * (4*pi/Area)^(3/2);

disp(['theoretical volume: ' num2str(Vol)]);
disp(['theoretical area: ' num2str(Area)]);
disp(['theoretical reduced volume: ' num2str(Red_Vol)]);

[x, y, z] = meshgrid(linspace(-250,250,64)); % simulation domain in nm
%F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;
F = sqrt(x.^2+y.^2+z.^2) - 215;
map = SD.SDF3(x,y,z,F);
map.reinitialization( map.F )
map.plotSurface(0,1,'g')

% calculate the Dirac delta function and the Heaviside function
Epsilon = map.GD3.Dx * 1.5;
Elt_Vol = map.GD3.Dx ^ 3;
In_vsc = map.F < -Epsilon;
Ot_vsc = map.F > Epsilon;
Ed_vsc = ~In_vsc & ~Ot_vsc;

Dirac_delta = zeros(map.GD3.Size);
Dirac_delta(Ed_vsc) =  map.Fg_1(Ed_vsc) .* (1 + cos(pi*map.F(Ed_vsc)/Epsilon)) / Epsilon / 2;

Heaviside = zeros(map.GD3.Size);
Heaviside(Ot_vsc) = 1;
Heaviside(Ed_vsc) = (1 + map.F(Ed_vsc)/ Epsilon + sin(pi*map.F(Ed_vsc)/Epsilon)/pi) / 2;

Cht_Min = 1 - Heaviside; % charateristic function of the interior region

Num_Vol = sum(Cht_Min(:)) .* Elt_Vol;
Num_Ara = sum(Dirac_delta(:)) .* Elt_Vol;
Num_Red_Vol = (3*Num_Vol/4/pi) * (4*pi/Num_Ara)^(3/2);

disp(['numerical volume: ' num2str(Num_Vol)]);
disp(['numerical area: ' num2str(Num_Ara)]);
disp(['numerical reduced volume: ' num2str(Num_Red_Vol)]);

disp(['relative error for volume: ' num2str(abs((Num_Vol-Vol)/Vol))]);
disp(['relative error for area: ' num2str(abs((Num_Ara-Area)/Area))]);

