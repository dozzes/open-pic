--[[
NAME
  em_field

FUNCTIONS
  uniform(grid, Ex, Ey, Ez, Bx, By, Bz)
  dipole_equator(grid, Bcloud, Bratio)
    
NOTES
    Package for EM field initialization
]]

local P = {}

if _REQUIREDNAME == nil then
    em_field = P
else
    _G[_REQUIREDNAME] = P
end


require("print_mpi")

function P.uniform(grid, Ex, Ey, Ez, Bx, By, Bz)
  
  for kx = 0, (grid:size_x() - 1) do
  for ky = 0, (grid:size_y() - 1) do
  for kz = 0, (grid:size_z() - 1) do
    grid:at(kx, ky, kz).E.x = Ex
    grid:at(kx, ky, kz).E.y = Ey
    grid:at(kx, ky, kz).E.z = Ez
    grid:at(kx, ky, kz).B.x = Bx
    grid:at(kx, ky, kz).B.y = By
    grid:at(kx, ky, kz).B.z = Bz
  end
  end
  end
end  -- set_uniform_field

-- Bz_ratio - Bmax/Bmin along z-axis
function P.dipole_equator(grid, Bcloud, Bratio)
  local B_max = 0.0
  local B_min = Bcloud
  
  local h = grid.step
  
  local cloud_x_node = math.floor(0.5 * grid:size_x())
  local cloud_y_node = math.floor(0.5 * grid:size_y())
  local cloud_z_node = math.floor(0.5 * grid:size_z())
                        
  local cloud_x = h * cloud_x_node
  local cloud_y = h * cloud_y_node
  local cloud_z = h * cloud_z_node
  
  local Lx = grid:length_x()
  local B_ratio_m3 =  math.pow(Bratio, 1.0/3.0)
  local p = -0.125*Bcloud*Lx*Lx*Lx*math.pow(((B_ratio_m3 + 1)/(B_ratio_m3 - 1)), 3.0)
    
  local dipole_x = -Lx /(B_ratio_m3 - 1.0)
    
  for kx = 0, (grid:size_x() - 1) do
  for ky = 0, (grid:size_y() - 1) do
  for kz = 0, (grid:size_z() - 1) do     
    
  	grid:at(kx, ky, kz).E.x = 0.0
    grid:at(kx, ky, kz).E.y = 0.0
    grid:at(kx, ky, kz).E.z = 0.0
    
    local x = kx*h - dipole_x
    local y = ky*h - cloud_y
    local z = kz*h - cloud_z
    
    local r2 = (x + 0.5*h) * (x + 0.5*h) + y*y + z*z
    local r5 = math.pow(r2, 2.5)
    local Bx = 3.0 *(x + 0.5*h) * z * p/r5
    grid:at(kx, ky, kz).B.x = Bx
 
    r2 = x*x + (y + 0.5*h) * (y + 0.5*h) + z*z
    r5 = math.pow(r2, 2.5)
    local By =  3.0*(y + 0.5*h) * z * p/r5
    grid:at(kx, ky, kz).B.y = By
    
    r2 = x*x + y*y + (z + 0.5*h)*(z + 0.5*h)
    r5 = math.pow(r2, 2.5)
    local Bz = (3.0 * (z + 0.5*h) * (z + 0.5*h) - r2)*p/r5
    grid:at(kx, ky, kz).B.z = Bz
    
    local B = math.sqrt(Bx*Bx + By*By + Bz*Bz)
    B_max = math.max(B, B_max)
  end
  end
  end
  return B_max
end  -- set_outer_magnetic_dipole_field
 