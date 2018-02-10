% test the mexReinitialization scheme

addpath(genpath('mexReinitialization'))
%path_test % test if  I can call function from another folder



	Kappa = 411;
	RV = 65; % reduced volume

	a = 215;
	b = 215;
	c = 62.3;

	xv = linspace(-250,250,64);
	yv = xv;
	%zv = xv(abs(xv)<100);
	zv = xv;

	[x, y, z] = meshgrid(xv, yv, zv); % simulation domain in nm

	F = sqrt(x.^2/a^2 + y.^2/b^2 + z.^2/c^2) - 1;

	map = SD.SDF3(x,y,z,F);
%	tic;map.reinitialization( map.F );toc

	%disp('calling mexReinitialization ...');

	% need to cast array into int32 corresponding to int in C
	%shift_mat = struct('soXo', int32(map.GD3.soXo), ...
	%				   'soxo', int32(map.GD3.soxo), ...
	%				   'sYoo', int32(map.GD3.sYoo), ...
	%				   'syoo', int32(map.GD3.syoo), ...
	%				   'sooZ', int32(map.GD3.sooZ), ...
	%				   'sooz', int32(map.GD3.sooz));

	% because c index start from 0, minus 1 is needed
	% index in c is a int32, conversion is needed
	shift_mat = struct('soXo', int32(map.GD3.soXo-1), ...
					   'soxo', int32(map.GD3.soxo-1), ...
					   'sYoo', int32(map.GD3.sYoo-1), ...
					   'syoo', int32(map.GD3.syoo-1), ...
					   'sooZ', int32(map.GD3.sooZ-1), ...
					   'sooz', int32(map.GD3.sooz-1));

	%tic
	%out=mexReinitialization(map.F, shift_mat,[map.GD3.Dx,map.GD3.Dy,map.GD3.Dz]);
	%toc

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 map.plotSurface(0,1,'g')

 loops = 0;
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

%
tic
 for ii = 1:loops-1
 %for ii = 1:0
	%tic;

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

	%disp('p1 time');
	%toc;
	

	if (mod(ii,SkipR)==0)
		cur_vol_BR = map.VolumeIntegral(1);
		cur_ara_BR = map.SurfaceIntegral(1);
		disp(['area error before reinitialization ', num2str(ii), ': ', num2str(cur_ara_BR/map.Suf_Area)]);
		disp(['volume error before reinitialization ', num2str(ii), ': ', num2str(cur_vol_BR/map.En_Volume)]);
	%tic;
		%map.reinitialization( reshape(F_new, map.GD3.Size) );
		map.F = mexReinitialization(map.F, shift_mat,[map.GD3.Dx,map.GD3.Dy,map.GD3.Dz]);
	%disp('reini time');
	%toc;
	end

	%tic;
	if (mod(ii,Skip)==0)
	%if false
	%	mov(count) = getframe(gcf);
		figure(ves)
		clf
		map.plotSurface(0,1,'g')
		time = num2str(ii*Dt);
		title([num2str(ii) ': ' num2str(ii*Dt)])
		text(map.GD3.xmin,map.GD3.ymax,(map.GD3.zmax+map.GD3.zmin)/2,['BR',num2str(ii),':',time])
		drawnow

		DistanceMap = map.F;
		%saveas(gcf, fullfile(Pic,[num2str(ii),'AA_BR','.png']))
		%saveas(gcf, fullfile(Pic,[num2str(ii),'AA_BR','.fig']))
		%save(fullfile(Mat,['DFV',num2str(ii),'AA_BR','.mat']),'DistanceMap')
	%	count = count + 1;
	end

	Eng_it = [Eng_it, cur_eng/Eng];
	vol_it = [vol_it, cur_vol/map.En_Volume-1]; 
	ara_it = [ara_it, cur_ara/map.Suf_Area-1];


	
	drawnow

	if (mod(ii,Skip)==0)
	%if false
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

		%saveas(gcf, fullfile(Pic,['eng','.png']))
	end

	disp(['bending energy ratio ', num2str(ii), ': ', num2str(c11/Eng)]);
	disp(['area error ', num2str(ii), ': ', num2str(cur_ara/map.Suf_Area)]);
	disp(['volume error ', num2str(ii), ': ', num2str(cur_vol/map.En_Volume)]);

	%keyboard

%	disp('p2 time');
%	toc
end

toc