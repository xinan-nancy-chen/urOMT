
function [smoothed_array] = affine_diffusion_3d(original_array,t_tot,dt,aff_flag,verbose)
% affine_diffusion_3d smooth an image using the affine-invariant
% mean_curvature flow. Runtime is a bit restrictive.
%
%   t_tot    = total evolution time, in "seconds"
%   dt       = numerical evolution parameter
%   aff_flag = 1 if affine-invariant flow, 0 for linear smoothing
%   verbose  = 1 if step counter and total runtime should be printed
%       (probably leave these flags both 1)
%
% Example usage:
% C = affine_diffusion_3d(A, .5, 0.02, 1, 1);

phi = double(original_array);    % cast to double to avoid rounding

%% Set up parameters
%dt = 0.01;
if ndims(phi) == 2 %#ok<ISMAT>
    error('Error: takes a 3D array.')
end
[Ny,Nx,Nz] = size(phi);


n_t = round(t_tot/dt);
%% Evolve the image
if(verbose)
    tic
end

for t = 1:n_t
    pX = zeros(size(phi)); pY = zeros(size(phi)); pZ = zeros(size(phi));
    pXX = zeros(size(phi)); pYY = zeros(size(phi)); pZZ = zeros(size(phi));
    if(aff_flag)
        pXY = zeros(size(phi)); pXZ = zeros(size(phi)); pYZ = zeros(size(phi));
    end
    
    % Compute derivatives for all voxels
    for j = 2:Nx-1    % x-dimension (columns)
        for k = 2:Ny-1    % y-dimension (rows)
            for l = 2:Nz-1    % z-dimension (pages) phi(k,j,l)
                pX(k,j,l) = 0.5*(phi(k,j+1,l)-phi(k,j-1,l));
                pY(k,j,l) = 0.5*(phi(k-1,j,l)-phi(k+1,j,l));
                pZ(k,j,l) = 0.5*(phi(k,j,l-1)-phi(k,j,l+1));
                pXX(k,j,l) = (phi(k,j+1,l)-2*phi(k,j,l)+phi(k,j-1,l));
                pYY(k,j,l) = (phi(k+1,j,l)-2*phi(k,j,l)+phi(k-1,j,l));
                pZZ(k,j,l) = (phi(k,j,l+1)-2*phi(k,j,l)+phi(k,j,l-1));
                if(aff_flag)
                    pXY(k,j,l) = 0.25*(-phi(k-1,j-1,l)+phi(k-1,j+1,l)...
                        + phi(k+1,j-1,l) - phi(k+1,j+1,l));
                    pXZ(k,j,l) = 0.25*(phi(k,j-1,l+1) + phi(k,j+1,l-1)...
                        - phi(k,j-1,l-1)-phi(k,j+1,l+1));
                    pYZ(k,j,l) = 0.25*(phi(k-1,j,l-1)+phi(k+1,j,l+1)...
                        -phi(k-1,j,l+1)-phi(k+1,j,l-1));
                end
            end
        end
    end
    
    if(aff_flag)
        % Compute curvature quantities
        meanCurvNum = pX.^2.*(pYY+pZZ) + pY.^2.*(pXX+pZZ) + pZ.^2.*(pXX+pYY) ...
            -2*(pX.*pY.*pXY + pX.*pZ.*pXZ + pY.*pZ.*pYZ);
        
        gausCurvNum = pX.^2.*(pYY.*pZZ-pYZ.^2) + pY.^2.*(pXX.*pZZ - pXZ.^2) ...
            + pZ.^2.*(pXX.*pYY - pXY.^2) + 2*pX.*pY.*(pXZ.*pYZ-pXY.*pZZ) ...
            + 2*pY.*pZ.*(pXY.*pXZ-pYZ.*pXX) + 2*pX.*pZ.*(pXY.*pYZ-pXZ.*pYY);
        
        % Update phi according to mean curvature flow
        %phi = phi + dt*sign(meanCurvNum).*nthroot(max(0,gausCurvNum),4);
        phi = max(zeros(size(phi)),phi + dt*sign(meanCurvNum).*nthroot(max(0,gausCurvNum),4));
    else % just do linear
        phi = phi + dt*(pXX+pYY+pZZ);
    end
    
    if(verbose)
        fprintf('Step %d out of %d.\n', t,n_t);
    end
end

if(verbose)
    rt = toc;
    fprintf('\n Evolution took %f seconds.\n', rt);
end
smoothed_array = phi;  % return

