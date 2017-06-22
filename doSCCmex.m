function doSCCmex(dist,debug)


if (exist('dist','var') ~= 1)
   dist = 0;
end
if (exist('debug','var') ~=1)
   debug = 0;
end

mex_files = {'InnerSCCmex.c' ...
   };
doGENmex(mex_files,dist,debug);
