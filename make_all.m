% Compiles MEX dependencies

% Build Graph-cuts solver
cd 3rd-party/max-flow
make
cd ../..

% Build Mex routines from the UGM toolbox
cd 3rd-party/ugm
mexAll
cd ../..

% Build DWT routines
cd 3rd-party/dde-prnu/Filter
compile
cd ../../..