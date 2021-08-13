%%%ExampleScript to calculate weighted control energy

T = 1;
rho = 1;

% Control nodes selection
n = 232;
xc = eye(n);

% load sum connectivity for integrated and segregated states;
integratedTable = readtable("C:\Users\Angelika\Dropbox\PhD\EXPERIMENTS\03_fMRI\ANALYSES\AverageStates\IntegratedSumConnectivity.csv");
segregatedTable = readtable("C:\Users\Angelika\Dropbox\PhD\EXPERIMENTS\03_fMRI\ANALYSES\AverageStates\SegregatedSumConnectivity.csv");

Network = load("C:\Users\Angelika\Dropbox\PhD\EXPERIMENTS\03_fMRI\ANALYSES\NBS\Labels\DMN_232.txt")

% % Nodes to be constrained
% S = zeros(n);
% for i = 1:n
%   if ismembertol(i, Network) == 1
%     S(i, i) = 1;
%   end
% end

% Constrain whole-brain
S = eye(n);

% loop through all subjects
all = [007, 010, 014, 015, 017, 028, 038, 039, 051, 056, 057, 068, 070, 073, 078, 087, 091, 092, 093, 099, 104, 108, 111, 117, 118, 119, 122, 124, 126, 128, 131, 132, 134, 137, 138, 142, 145, 146, 149, 150, 152, 156, 157, 158, 159, 160, 162, 163, 167, 169, 172, 175, 178, 179, 184, 185, 186, 187, 188, 189, 190, 193, 194, 195, 197, 198, 199, 201, 208, 209, 210, 214, 215, 216, 220, 222, 223, 227, 228, 229, 232, 233, 234, 235, 236, 238, 239, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 256, 257, 264, 268, 272, 274, 275, 276, 277, 280, 283, 284, 309, 311, 377, 404, 504, 601 ];


for s = 1:numel(all)
	name = sprintf('VIPD_%03d', all(s));
	ConnPath = fullfile("C:\Users\Angelika\Dropbox\PhD\EXPERIMENTS\03_fMRI\ANALYSES\Structure-Function\CONNECTOMES", sprintf('VIPD_%03d_shaefer232_connectome.csv', all(s)));
	x0 = integratedTable.(name);
	xf = integratedTable.(name);
	outFile = fullfile("C:\Users\Angelika\Dropbox\PhD\EXPERIMENTS\03_fMRI\ANALYSES\NetworkControl\", 'IntegratedPersistence', name);
    EnergyCal_Function(ConnPath, T, xc, x0, xf, S, rho, outFile);
end
