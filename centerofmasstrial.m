% Center of mass in one dimension
% Modifiled by Song, 2024.11.16, Chongqing University

clear; % Delete all the varibles
ruler.length = 30;   % cm
ruler.weight = 25;   % gram

mesh.dx = 0.1;
mesh.xs = 0:mesh.dx:ruler.length;   % coordinates in cm

ruler.density = ruler.weight/numel(mesh.xs); % grams / cell


masses = ones(size(mesh.xs))*ruler.density;
M = sum(masses);

COM = centerOfMass(masses,mesh);

mass = 100;    % gram
position = 3;  % cm
masses = addMass(mass,position,masses,mesh);
COM = centerOfMass(masses, mesh); 
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Center of mass in two dimensions

plate.length = 14;   % cm
plate.width = 7;     % cm
plate.weight = 77;   % gram

clear mesh

mesh.dx = 0.1;
mesh.dy = 0.1;
mesh.xs = 0:mesh.dx:plate.length;
mesh.ys = 0:mesh.dy:plate.width;

[mesh.xgrid, mesh.ygrid] = meshgrid(mesh.xs,mesh.ys); %The grid vector mesh.xgrid is replicated numel(mesh.xs) times to form the columns of X. Same applies to Y

plate.density = plate.weight/numel(mesh.xgrid); 
masses_2D = plate.density * ones(size(mesh.xgrid));

M_2Dplate = matrixsum(masses_2D);

COM_2D = COMof2D(masses_2D, mesh);

COM_2D

weight_2D = 100; %gram
position_2D = [0.3,0.3]; %cm

masses_2D = Addmass_2D(weight_2D, position_2D, masses_2D, mesh);

new_weight = matrixsum(masses_2D);

 
COM2D_update_x = matrixsum(masses_2D.*mesh.xgrid)/new_weight;
COM2D_update_y = matrixsum(masses_2D.*mesh.ygrid)/new_weight;

COM2D_update = [COM2D_update_x, COM2D_update_y];

COM2D_update
% Compare
COM_2D_compare = COMof2D(masses_2D, mesh);

% plate
% Next we create a matrix that's the same size as the mesh, where the value at each location is the density.
% masses = ones(size(mesh.xgrid)) * plate.density
% The matrixSum function still works (in fact, it works with any number of dimensions).  We can use it to confirm that the total mass is correct.
% M = matrixSum(masses)
% Computing the center of mass in two dimensions takes a little more work.  See centerOfMass2 below and make sure you understand it.
% COM = centerOfMass2(masses,mesh)
% Does the result make sense?
% Now let's put a 100 gram mass at the position x=0.3 cm, y=0.3 cm, measured from a corner of the plate.  To represent a point in 2D, I could use a vector like this:
% position = [0.3, 0.3];     % cm
% Or a structure with fields x and y, as shown below.
% mass = 100;   % gram
% 
% clear position
% position.x = 0.3;   % cm
% position.y = 0.3;   % cm
% To add the mass to the plate, we need a new version of addMass2.
% Exercise 2: Fill in the body of addMass2, below, so it does what it is supposed to.
% masses = addMass2(mass,position,masses,mesh)
% Once you have addMass2 working, compute the COM again and see if it shifted as expected.  Check the result by hand.
% M = matrixSum(masses)
% COM = centerOfMass2(masses,mesh)
% 
% Exercise 3:  Check your notes from the last class to see how you arranged weights on your plate.  Modify the code above to represent the scenario you created.  Adjust the parameters of the plate if nessesary.  Run the code and see if the answer is consistent with what you calculated yesterday.
% 
% 
% Center of mass in three dimensions
% This is an optional exercise if you finish the previous sections and have some additional time.  It is a good opportunity to practice the MATLAB features we have learned so far and continue to develop your progamming skills.
% Exercise 4: Now that we have two dimensions, three dimensions is easy!
% Create a structure called cube that contains length = 10 cm, width = 10 cm, height = 10 cm, and weight = 1000 grams.
% Create a mesh with dx, dy, and dz = 0.2.
% Compute cube.density and create a matrix called masses with the same size as the mesh.  Use matrixSum to confirm that the total mass is right.
% Write a function called centerOfMass3 that computes the center of mass in three dimensions.  Test your function and confirm that the center of mass is where you expect it to be.
% Suppose we replace one cell of the cube, in one corner, with a super dense material, so that the mass of that cell is 1000 grams.  Where do you expect the new center of mass to be?  Compute it and see if you got it right.
% 
% 
% Function definitions
% function M = matrixSum(masses)
%     % matrixSum: returns total of all elements in the matrix
%     
%     % normally sum(m) computes the sums of the columns
%     % selecting m(:) flattens the matrix and computes the sum of all elements
%     % see https://stackoverflow.com/questions/1721987/what-are-the-ways-to-sum-matrix-elements-in-matlab
%     M = sum(masses(:));
% end
% 
% function COM = centerOfMass(masses,mesh)
%     % centerOfMass: computes center of mass in 1D
%     % masses: matrix of masses
%     % mesh: structure containing xs
%     % returns: scalar
%     M = matrixSum(masses);
%     COM = matrixSum(masses .* mesh.xs) / M;
% end
% 
% function masses = addMass(mass,position,masses,mesh)
%     % addMass: adds mass at a given position in 1D
%     % mass: scalar, in grams
%     % position: scalar, in cm
%     % masses: matrix of masses
%     % mesh: structure containing xs
%     % returns: the updated matrix of masses
%     
%     % find the index of the location in the mesh closest to position
%     index = closestIndex(position,mesh.xs);
%     
%     % update the matrix
%     masses(index) = masses(index) + mass;
% end
% 
% function index = closestIndex(coord,grid)
%     % closestIndex: finds the index of the grid point closest to coord
%     % coord: scalar coordinate
%     % grid: matrix of coordinates
%     % returns: integer index
%     
%     % see: https://www.mathworks.com/matlabcentral/answers/152301-find-closest-value-in-array
%     [c index] = min(abs(coord-grid));
% end
% 
% function COM = centerOfMass2(masses,mesh)
%     % centerOfMass2: computes center of mass in 2D
%     % masses: matrix of masses
%     % mesh: structure containing xgrid and xgrid
%     % returns: Vector [xcom,ycom]
%     
%     M = matrixSum(masses);
%     xcom = matrixSum(masses .* mesh.xgrid) / M;
%     ycom = matrixSum(masses .* mesh.ygrid) / M;
%     COM = [xcom,ycom];
% end
% 
% function masses = addMass2(mass,position,masses,mesh) 
%     % addMass2: adds mass at a given position in 2D
%     % mass: scalar, in grams
%     % position: scalar, in cm
%     % masses: matrix of masses
%     % mesh: structure containing xs
%     % returns: the updated matrix of masses
%     
%     % TODO: FIX THIS SO IT PUTS THE MASS IN THE RIGHT PLACE
%     masses(1,1) = masses(1,1) + mass;
% end


function M = matrixsum(masses)

M = sum(masses(:));

end



function COM = centerOfMass(masses,mesh)
    M = sum(masses);
    COM = sum(masses .* mesh.xs) / M;
end


function masses = addMass(mass,position,masses,mesh)
    % addMass: adds mass at a given position in 1D
    % mass: scalar, in grams
    % position: scalar, in cm
    % masses: matrix of masses
    % mesh: structure containing xs
    % returns: the updated matrix of masses
    
    % find the index of the location in the mesh closest to position
    index = closestIndex(position,mesh.xs);
    
    % update the matrix
    masses(index) = masses(index) + mass;
end

function index = closestIndex(coord,grid)
    [c index] = min(abs(coord-grid));
end


function COM2D = COMof2D(masses_2D, mesh) 

x = matrixsum(masses_2D.*mesh.xgrid)/matrixsum(masses_2D);
y = matrixsum(masses_2D.*mesh.ygrid)/matrixsum(masses_2D);

COM2D = [x, y];

end

function masses = Addmass_2D(weight, position, masses, mesh)

index_x = closestIndex(position(1,1), mesh.xgrid(1,:));
 
index_y = closestIndex(position(1,2), mesh.ygrid(:,1));
 
masses(index_x, index_y) = masses (index_x, index_y) + weight;

end

