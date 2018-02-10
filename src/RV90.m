% before running a test, gather necessary information about this test
% so that the result can be recovered easily

% time of excution
	Time_Now = datetime('now');
	FormatOut = 'yy_mm_dd_HH_MM_SS_';
	Time_Begin = datestr(Time_Now,FormatOut);

% system info
	if ismac
		PlatForm = 'Mac_';
		PaFold = '..';
	elseif isunix
		PlatForm = 'unix_';
		PaFold = '/extra/tiankuizhang';
	elseif ispc
		PlatForm = 'Windows_';
		PaFold = '..';
	else
		PlatForm = 'unknown_';
	end

% make a directory tagged with time of excution outsider src folder
	% path of this folder. we randomized the last 3 digits to avoid name clashes
	Result_Folder = fullfile(PaFold,'exc',[PlatForm Time_Begin num2str(floor(rand()*1000))]); 
	if ~exist(fullfile('..','exc')) % if ../exc does not exist
		mkdir(fullfile('..','exc')) % make this folder
	end
	while exist(Result_Folder) % if this folder name alread exists -- which will probably never happen
		Result_Folder = fullfile('..','exc',[PlatForm Time_Begin num2str(floor(rand()*1000))]); % rename it
	end
	mkdir(Result_Folder) % make a folder

% folder to hold png and mat files
	Pic = fullfile(Result_Folder, 'pic');
	Mat = fullfile(Result_Folder, 'mat');

	mkdir(Pic)
	mkdir(Mat)
	
% now copy the src code the Result_Folder
	Current_src = fullfile(Result_Folder,'src'); % folder to hold current src files
	mkdir(Current_src)
	copyfile(fullfile('..','src'),Current_src);

% now open a text file and write to it comments and testing info
	Test_info = fopen(fullfile(Result_Folder,'test_info'), 'w');
	fprintf(Test_info, 'test start time : \t');
	fprintf(Test_info, [datestr(Time_Now, 'yy/mm/dd HH:MM:SS'),'\n']);
	fprintf(Test_info, 'test file: VD \n');
	fprintf(Test_info, 'test info: try to accerlerate convergence of iters \n');
	fclose(Test_info);

% diary the command window
	diary(fullfile(Result_Folder,'log_command_window'))
	diary on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Kappa = 411;

RV = 90; % reduced volume

switch RV
	case 90

		a = 210;
		b = 210;
		c = 114;

		xv = linspace(-250,250,64);
		yv = xv;
		zv = xv(abs(xv)<150);

		[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm


		F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	case 80

		a = 200;
		b = 150;
		c = 71;

		xv = linspace(-250,250,64);
		yv = xv(abs(xv)<200);
		zv = xv(abs(xv)<100);

		[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm


		F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	case 70

		a = 215;
		b = 215;
		c = 70;

		xv = linspace(-250,250,64);
		yv = xv;
		zv = xv(abs(xv)<100);

		[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm


		F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	case 65

		a = 215;
		b = 215;
		c = 62.3;

		xv = linspace(-250,250,64);
		yv = xv;
		zv = xv(abs(xv)<100);

		[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm


		F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	case 60

		a = 215;
		b = 215;
		c = 55.5;

		xv = linspace(-250,250,64);
		yv = xv;
		zv = xv(abs(xv)<100);

		[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm


		F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	case 55

		a = 215;
		b = 215;
		c = 49.3;

		xv = linspace(-250,250,64);
		yv = xv;
		zv = xv(abs(xv)<100);

		[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm


		F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	case 50

		a = 215;
		b = 215;
		c = 43.6;

		xv = linspace(-250,250,64);
		yv = xv;
		zv = xv(abs(xv)<100);

		[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm


		F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	otherwise
		disp('unknown RV!');
end

% vesicle dynamics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%F = sqrt(x.^2+y.^2+z.^2) - 215;
map = SD.SDF3(x,y,z,F);
map.reinitialization( map.F )
%map.plotSurface(0,1,'g')

%keyboard

% save grids
GridX = map.GD3.X;
GridY = map.GD3.Y;
GridZ = map.GD3.Z;
save(fullfile(Result_Folder,'Grid.mat'),'GridX','GridY','GridZ');

loops = 2000;
Skip = 20;
SkipR = 1;
Dt = 2 * map.GD3.Dx ^ 4 / Kappa; % it is more like Dt*zeta(the drag coefficient)

Eng = map.SurfaceIntegral(map.SC.^2);

Eng_it = [];
vol_it = [];
ara_it = [];

%keyboard

ves = figure('Name', 'Vesicle Shape');
eng = figure('Name', 'energy/volume/area vs time');

for ii = 1:loops-1

	cur_eng = map.SurfaceIntegral(map.SC.^2);
	cur_vol = map.VolumeIntegral(1);
	cur_ara = map.SurfaceIntegral(1);

	err_ara = (map.Suf_Area - cur_ara) / Dt;
	err_vol = (map.En_Volume - cur_vol) / Dt;

	% calculate normal velocity due to bending
	co = 0.5 * (map.SC.^2 - 4*map.GC); % will be reused
	Vb = Kappa * (map.SL(map.SC) + map.SC .* co);
	% rate of chage change of area and volume due to Vb minus error
	DA = map.SurfaceIntegral(map.SC .* Vb) - err_ara;
	DV = map.SurfaceIntegral(Vb) - err_vol;
	%DA = (cur_ara - map.Suf_Area) / Dt;
	%DV = (cur_vol - map.En_Volume) / Dt;
	s_mat = [DA; DV];
	% coefficient matrix from lagrange multiplier
	%c11 = map.SurfaceIntegral(map.SC.^2);
	c11 = cur_eng;
	c12 = - map.SurfaceIntegral(map.SC);
	c21 = - c12;
	%c22 = - map.Suf_Area;
	c22 = - cur_ara;
	c_mat = [c11, c12; c21, c22]; 
	% calculate surface tension and pressure
	lag_mul = c_mat\s_mat; % lagrange multiplier
	Tension = lag_mul(1);
	Pressure = lag_mul(2);

	%% construct the force density operator
	% an operator appearing in surface laplacian operator and curvature operator
	LOP = map.GD3.Lxx + map.GD3.Lyy + map.GD3.Lzz ...
				- ( map.GD3.SparseDiag(map.Nx .* map.Nx) * map.GD3.Lxx + ...
					map.GD3.SparseDiag(map.Ny .* map.Ny) * map.GD3.Lyy + ...
					map.GD3.SparseDiag(map.Nz .* map.Nz) * map.GD3.Lzz + ...
					map.GD3.SparseDiag(map.Nx .* map.Ny) * map.GD3.Lxy * 2 + ...
					map.GD3.SparseDiag(map.Ny .* map.Nz) * map.GD3.Lyz * 2 + ...
					map.GD3.SparseDiag(map.Nz .* map.Nx) * map.GD3.Lzx * 2  ); 
	% curvature operator 	
	CvL = map.GD3.SparseDiag(1./map.Fg) * LOP;
	% normal derivative
	ND = map.GD3.SparseDiag(map.Nx) * map.GD3.Lx + ...
		 map.GD3.SparseDiag(map.Ny) * map.GD3.Ly + ...
		 map.GD3.SparseDiag(map.Nz) * map.GD3.Lz ;
	% surface laplacian operator
	SL = LOP - map.GD3.SparseDiag(map.SC) * ND;
	% normal velocity operator without drag coefficient and Pressure
	NV = Kappa * ( ( SL + map.GD3.SparseDiag(co) ) * CvL ) - Tension * CvL;
	%NV = Kappa * ( ( SL + map.GD3.SparseDiag(co) ) * CvL );
	%NV = Kappa * SL  * CvL;

	% the advection term without pressure
	A = map.GD3.SparseDiag(map.Fg_1) * NV * Dt / 2;
	B = map.GD3.Idt + A;
	C = map.GD3.Idt - A;

	F_old = map.F;
	S = C * F_old(:) - Dt * Pressure * map.Fg_1(:);
	%S = C * F_old(:);

	[L,U]=ilu(B,struct('type','nofill','milu','row'));
	%F_new = gmres(B, S, [], 1e-12, 300, L, U);
	F_new = bicgstab(B, S, 1e-12, 300, L, U);
	%[F_new, flag] = bicgstab(B, S, 1e-12, 300, L, U);
	%if flag
	%	keyboard
	%end
	map.F = reshape(F_new, map.GD3.Size);

	

	if (mod(ii,SkipR)==0)
		cur_vol_BR = map.VolumeIntegral(1);
		cur_ara_BR = map.SurfaceIntegral(1);
		disp(['area error before reinitialization ', num2str(ii), ': ', num2str(cur_ara_BR/map.Suf_Area)]);
		disp(['volume error before reinitialization ', num2str(ii), ': ', num2str(cur_vol_BR/map.En_Volume)]);
		map.reinitialization( reshape(F_new, map.GD3.Size) );
	end

	if (mod(ii,Skip)==0)
	%	mov(count) = getframe(gcf);
		figure(ves)
		clf
		map.plotSurface(0,1,'g')
		time = num2str(ii*Dt);
		title([num2str(ii) ': ' num2str(ii*Dt)])
		text(map.GD3.xmin,map.GD3.ymax,(map.GD3.zmax+map.GD3.zmin)/2,['BR',num2str(ii),':',time])
		drawnow

		DistanceMap = map.F;
		saveas(gcf, fullfile(Pic,[num2str(ii),'AA_BR','.png']))
		%saveas(gcf, fullfile(Pic,[num2str(ii),'AA_BR','.fig']))
		save(fullfile(Mat,['DFV',num2str(ii),'AA_BR','.mat']),'DistanceMap')
	%	count = count + 1;
	end

	Eng_it = [Eng_it, cur_eng/Eng];
	vol_it = [vol_it, cur_vol/map.En_Volume-1]; 
	ara_it = [ara_it, cur_ara/map.Suf_Area-1];


	
	drawnow

	if (mod(ii,Skip)==0)
		figure(eng)
		clf

		yyaxis left
		plot(1:ii, Eng_it)
		ylabel('energy measusred in initial energy')
		xlabel('iteration')

		yyaxis right
		plot(1:ii, vol_it)
		hold on
		plot(1:ii, ara_it)
		hold off
		ylabel('deviation from initial area and volume')

		legend('energy','volume error','area error')

		drawnow

		saveas(gcf, fullfile(Pic,['eng','.png']))
	end



	disp(['bending energy ratio ', num2str(ii), ': ', num2str(c11/Eng)]);
	disp(['area error ', num2str(ii), ': ', num2str(cur_ara/map.Suf_Area)]);
	disp(['volume error ', num2str(ii), ': ', num2str(cur_vol/map.En_Volume)]);

	%keyboard
end

% write test end time

	Test_info = fopen(fullfile(Result_Folder,'test_end_info'), 'w');
	fprintf(Test_info, 'test end time : \t');
	fprintf(Test_info, [datestr(datetime('now'), 'yy/mm/dd HH:MM:SS'),'\n']);
	fclose(Test_info);
	diary off
