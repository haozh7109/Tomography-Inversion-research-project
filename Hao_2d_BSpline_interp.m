function [data_intp,data_x1y0,data_x0y1,data_x2y0,data_x0v2,data_x1v1] = Hao_2d_BSpline_interp(data_input,grid_x,grid_y,interp_loc)
%Job used to 2D interpolaiton based on 2d bi-cubic spline interpolation
%
% data_input:  input data for the interpolation use.
% grid_x:      X coordinate of the input data grid.
% grid_y:      y coordinate of the input data grid.
% interp_loc:  X,y coordinate of the location to be interpolated with format: [x , y].

%Coded by Hao, Nov 10. 2016
%% set default parameters used for function demo
if(nargin<1)
    
    %===== test on gaussian model ==========
    load('test_gaussian.mat')
    grid_x = x;
    grid_y = y;
    data_input = Vint;
    
    dx = grid_x(2)- grid_x(1);
    dy = grid_y(2)- grid_y(1);
    
    %define a interpolation grid
    n_grid_x = min(grid_x)+dx*0:dx/2:max(grid_x)-dx*0;
    n_grid_y = min(grid_y)+dy*0:dy/2:max(grid_y)-dy*0;
    
end

%% start the interpolation test

%derive the current coorinates
cod = interp_loc;

%% basic parameter setting
x  = cod(1);
y  = cod(2);
dx = grid_x(2) - grid_x(1);
dy = grid_y(2) - grid_y(1);


%% locate the current calculation point
ix0 = floor((x -grid_x(1))/dx)+1;
iy0 = floor((y -grid_y(1))/dy)+1;

%% calculate the dimontionless vector u(x) and v(y)
u = (x - grid_x(ix0))/dx;
v = (y - grid_y(iy0))/dy;

%% generate the cubic basis function for u and v
% ---derive the 0 order derivatives
bu_0_m1 = 1/6 * (1-u)^3;
bu_0_m0 = 1/6 * (3*u^3  - 6*u^2 + 4 );
bu_0_p1 = 1/6 * (-3*u^3 + 3*u^2 + 3*u + 1);
bu_0_p2 = 1/6 * (u)^3;
bu_0 = [bu_0_m1,bu_0_m0,bu_0_p1,bu_0_p2];

bv_0_m1 = 1/6 * (1-v)^3;
bv_0_m0 = 1/6 * (3*v^3  - 6*v^2 + 4 );
bv_0_p1 = 1/6 * (-3*v^3 + 3*v^2 + 3*v + 1);
bv_0_p2 = 1/6 * (v)^3;
bv_0 = [bv_0_m1,bv_0_m0,bv_0_p1,bv_0_p2];

% ---derive the 1 order derivatives
bu_1_m1 = 1/2 * (-u^2   + 2*u - 1);
bu_1_m0 = 1/2 * (3*u^2  - 4*u );
bu_1_p1 = 1/2 * (-3*u^2 + 2*u + 1);
bu_1_p2 = 1/2 * (u)^2;
bu_1 = [bu_1_m1,bu_1_m0,bu_1_p1,bu_1_p2];

bv_1_m1 = 1/2 * (-v^2   + 2*v - 1);
bv_1_m0 = 1/2 * (3*v^2  - 4*v );
bv_1_p1 = 1/2 * (-3*v^2 + 2*v + 1);
bv_1_p2 = 1/2 * (v)^2;
bv_1 = [bv_1_m1,bv_1_m0,bv_1_p1,bv_1_p2];

% ---derive the 2 order derivatives
bu_2_m1 =   -u  + 1;
bu_2_m0 =  3*u  - 2;
bu_2_p1 = -3*u  + 1;
bu_2_p2 =    u;
bu_2 = [bu_2_m1,bu_2_m0,bu_2_p1,bu_2_p2];

bv_2_m1 =   -v  + 1;
bv_2_m0 =  3*v  - 2;
bv_2_p1 = -3*v  + 1;
bv_2_p2 =    v;
bv_2 = [bv_2_m1,bv_2_m0,bv_2_p1,bv_2_p2];

%% extract the data along the grids
%pad the input data (padding 2 extra line on all boundaries)
[r_org,c_org]     = size(data_input);
data_pad          = zeros((r_org+4),(c_org+4));
data_pad(3:r_org+2,3:c_org+2) = data_input;
data_pad(1,:)     = data_pad(3,:);
data_pad(2,:)     = data_pad(3,:);
data_pad(end,:)   = data_pad(end-2,:);
data_pad(end-1,:) = data_pad(end-2,:);
data_pad(:,1)     = data_pad(:,3);
data_pad(:,2)     = data_pad(:,3);
data_pad(:,end)   = data_pad(:,end-2);
data_pad(:,end-1) = data_pad(:,end-2);

%extract the data from grid, according to current point location
% M_uv = zeros(4,4);
% M_uv = data_input(iy0-1:1:iy0+2,ix0-1:1:ix0+2);
M_uv = data_pad(iy0+1:1:iy0+4,ix0+1:1:ix0+4);

%% do the interpolation
sum_u0_v0 = 0;
sum_u1_v0 = 0;
sum_u0_v1 = 0;
sum_u2_v0 = 0;
sum_u0_v2 = 0;
sum_u1_v1 = 0;

for iv= 1:4
    
    %set the temporay summation.
    sum_u0   = 0;
    sum_u1   = 0;
    sum_u0v  = 0;
    sum_u2   = 0;
    sum_u0v2 = 0;
    sum_u1v1 = 0;
    
    %start summation from u(x) direction.
    for iu= 1:4
        %(1)derive the interpolated value
        sum_u0   =  sum_u0   + bu_0(iu) * M_uv(iv,iu);
        %(2)derive the first order derivative with respect to u(x)
        sum_u1   =  sum_u1   + bu_1(iu) * M_uv(iv,iu);
        %(3)derive the first order derivative with respect to v(y)
        sum_u0v  =  sum_u0v  + bu_0(iu) * M_uv(iv,iu);
        %(4)derive the second order derivative with respect to u(x)
        sum_u2   =  sum_u2   + bu_2(iu) * M_uv(iv,iu);
        %(5)derive the second order derivative with respect to v(y)
        sum_u0v2 =  sum_u0v2 + bu_0(iu) * M_uv(iv,iu);
        %(6)derive the mix partial derivative
        sum_u1v1 =  sum_u1v1 + bu_1(iu) * M_uv(iv,iu);
    end
    
    %full summation
    sum_u0_v0 = sum_u0_v0 + sum_u0   * bv_0(iv) *  1;
    sum_u1_v0 = sum_u1_v0 + sum_u1   * bv_0(iv) * (1/dx);
    sum_u0_v1 = sum_u0_v1 + sum_u0v  * bv_1(iv) * (1/dy);
    sum_u2_v0 = sum_u2_v0 + sum_u2   * bv_0(iv) * (1/dx)^2;
    sum_u0_v2 = sum_u0_v2 + sum_u0v2 * bv_2(iv) * (1/dy)^2;
    sum_u1_v1 = sum_u1_v1 + sum_u1v1 * bv_1(iv) * (1/dx)* (1/dy);
    
end

% saving the interpolated result and corresponding derivatives.
data_intp=sum_u0_v0;
data_x1y0=sum_u1_v0;
data_x0y1=sum_u0_v1;
data_x2y0=sum_u2_v0;
data_x0v2=sum_u0_v2;
data_x1v1=sum_u1_v1;

end


