% script to compile 

% -v : display detailed build and troubleshooting information
% -output nameofmex : control name of the mex file (for instance: mex -output mymex myfunc.c)
% -outdir path_to_dir : control the path to output dir
% to add -Wall flag, use: mex -v COMPFLAGS-'$cCOMPFLAGS -Wall' yprime.c
% -n : preview without executing 
% -c : compiler object file only
% -O : optimizes the object code

f1 = fullfile('mexReinitialization','mexReinitialization.c');
f2 = fullfile('mexReinitialization','Reinitialization.c');
f3 = fullfile('mexReinitialization','reinitialization_step.c');

%mex(f1,f2)

%compiler_option = '/GL /fp:fast -std=c99';
%compiler_option = '-std=c99';
%mex(f1,f2,['COMPFLAGS="$COMPFLAGS ' compiler_option '"']);

if ismac
	mex(f1,f2,f3)
elseif isunix
	mex(f1,f2,f3)
elseif ispc
	%mex('CFLAGS="\$CFLAGS -std=c99"',f1,f2)
	%mex -v CXXFLAGS='$CXXFLAGS -std=gnu99' ...
	%	mexReinitialization/mexReinitialization.c ...
	%	mexReinitialization/Reinitialization.c
	mex(f1,f2,f3,['CXXFLAGS="$CXXFLAGS -Wall"'])
end
	
	
