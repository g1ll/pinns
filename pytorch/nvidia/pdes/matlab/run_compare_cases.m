format long;
root = fileparts(mfilename('fullpath'));
result_root = fullfile(root, 'results');
if ~exist(result_root, 'dir')
	mkdir(result_root);
end
case1_dir = fullfile(result_root, 'case1');
case2_dir = fullfile(result_root, 'case2');
if ~exist(case1_dir, 'dir')
	mkdir(case1_dir);
end
if ~exist(case2_dir, 'dir')
	mkdir(case2_dir);
end
t1 = dof2y(1,20,0.14,14,0.45,0,0,0.02,0.02,0.02,1,50,fullfile(case1_dir,'case1'),1,0,'Caso 1 - Forma em T');
fprintf('CASE1_TMAX=%.15f\n', t1);
t2 = dof2y(1,10,0.07,0.07,0.5,13.6875,13.6875,0.015,0.015,0.04,1,50,fullfile(case2_dir,'case2'),1,0,'Caso 2 - Forma otima');
fprintf('CASE2_TMAX=%.15f\n', t2);