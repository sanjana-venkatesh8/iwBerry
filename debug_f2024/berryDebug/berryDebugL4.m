%% decode location and orientation of berry L4's
% MJB: index_L4 = orientation + num_orient * x_loc + num_orient * num_x * y_loc

berryL4inputs = readmatrix("Jens_input_L4_save0.txt");

L4Orients = mod(berryL4inputs, 4); % number of orientations = 4;
inter = (berryL4inputs - L4Orients) / 4; % xLoc + nX*yLoc;
L4XLoc = mod(inter, 8); % nX = 8
L4YLoc = (inter - L4XLoc) / 8;

% add 1 to all bc of off by 1 indexing
L4Orients = L4Orients + 1;
L4XLoc = L4XLoc + 1;
L4YLoc = L4YLoc + 1;

% change to my encoding for L4 inputs
iOrientedL4 = (L4Orients - 1) * 8 * 8;
% set the corresponding L4 neurons to ON
berryToVenkateshL4Indices = iOrientedL4 + ((L4YLoc - 1) * 8 + L4XLoc);

writematrix(berryToVenkateshL4Indices, "b2VL4Indices.txt");