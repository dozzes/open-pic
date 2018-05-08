--[[
NAME
  bound_cond

FUNCTIONS
   nonperturbed_NP(grid)
   nonperturbed_UP(grid)
   nonperturbed_UE(grid)
   nonperturbed_MF(grid)
   nonperturbed_EF(grid)
    
NOTES
    Package for boundary confitions
]]

local P = {}

if _REQUIREDNAME == nil then
    bound_cond = P
else
    _G[_REQUIREDNAME] = P
end

function P.nonperturbed_NP(grid)
  local to_x = grid:size_x() - 1
  local to_y = grid:size_y() - 1
  local to_z = grid:size_z() - 1
  
  local NP = 0.0
  
  -- Z 
  for kx = 0, to_x do
  for ky = 0, to_y do
    NP = grid:at(kx, ky, (to_z - 3)).NP                                         
    grid:at(kx, ky, (to_z - 2)).NP, grid:at(kx, ky, (to_z - 1)).NP, grid:at(kx,ky,to_z).NP = NP, NP, NP
       
    NP = grid:at(kx, ky, 2).NP; 
    grid:at(kx, ky, 1).NP, grid:at(kx, ky, 0).NP = NP, NP
  end
  end
  
  -- Y
  for kx = 0, to_x do
  for kz = 0, to_z do
    NP = grid:at(kx, (to_y - 3), kz).NP
    grid:at(kx, (to_y - 2), kz).NP, grid:at(kx, (to_y - 1), kz).NP, grid:at(kx, to_y, kz).NP = NP, NP, NP
       
    NP = grid:at(kx, 2, kz).NP;
    grid:at(kx, 1, kz).NP, grid:at(kx, 0, kz).NP = NP, NP
  end
  end
  
  --X
  for ky = 0, to_y do
  for kz = 0, to_z do
    NP = grid:at((to_x - 3), ky, kz).NP 
    grid:at((to_x - 2), ky, kz).NP, grid:at((to_x - 1), ky, kz).NP, grid:at(to_x, ky, kz).NP = NP, NP, NP
       
    NP = grid:at(2, ky, kz).NP
    grid:at(1, ky, kz).NP, grid:at(0, ky, kz).NP = NP, NP
  end
  end
end

function P.nonperturbed_UP(grid)
  local to_x = grid:size_x() - 1
  local to_y = grid:size_y() - 1
  local to_z = grid:size_z() - 1
  
  local UPx = 0.0
  local UPy = 0.0
  local UPz = 0.0
  
  -- Z 
  for kx = 0, to_x do
  for ky = 0, to_y do
  	UPx = grid:at(kx, ky, (to_z - 3)).UP.x
  	grid:at(kx, ky, (to_z - 2)).UPx, grid:at(kx, ky, (to_z - 1)).UP.x, grid:at(kx, ky, to_z).UP.x = UPx, UPx, UPx
    
    UPy = grid:at(kx, ky, (to_z - 3)).UP.y
    grid:at(kx, ky, (to_z - 2)).UP.y, grid:at(kx, ky, (to_z - 1)).UP.y, grid:at(kx, ky, to_z).UP.y = UPy, UPy, UPy
    
    UPz = grid:at(kx, ky, (to_z - 3)).UP.z
    grid:at(kx, ky, (to_z - 2)).UP.z, grid:at(kx, ky, (to_z - 1)).UP.z, grid:at(kx, ky, to_z).UP.z = UPz, UPz, UPz

    UPx = grid:at(kx, ky, 2).UP.x
    grid:at(kx, ky, 1).UP.x, grid:at(kx, ky, 0).UP.x = UPx, UPx
    
    UPy = grid:at(kx, ky, 2).UP.y
    grid:at(kx, ky, 1).UP.y, grid:at(kx, ky, 0).UP.y = UPy, UPy
    
    UPz = grid:at(kx, ky, 2).UP.z
    grid:at(kx, ky, 1).UP.z, grid:at(kx, ky, 0).UP.z = UPz, UPz
  end
  end
  
  -- Y
  for kx = 0, to_x do
  for kz = 0, to_z do
    UPx = grid:at(kx, (to_y - 3), kz).UP.x
    grid:at(kx, (to_y - 2), kz).UP.x, grid:at(kx, (to_y - 1), kz).UP.x, grid:at(kx, to_y, kz).UP.x = UPx, UPx, UPx
    
    UPy = grid:at(kx, (to_y - 3), kz).UP.y
    grid:at(kx, (to_y - 2), kz).UP.y, grid:at(kx, (to_y - 1), kz).UP.y, grid:at(kx, to_y, kz).UP.y = UPy, UPy, UPy
    
    UPz = grid:at(kx, (to_y - 3), kz).UP.z
    grid:at(kx, (to_y - 2), kz).UP.z, grid:at(kx, (to_y - 1), kz).UP.z, grid:at(kx, to_y, kz).UP.z = UPz, UPz, UPz
                   
    UPx = grid:at(kx, 2, kz).UP.x
    grid:at(kx, 1, kz).UPx, grid:at(kx, 0, kz).UP.x = UPx, UPx
    
    UPy = grid:at(kx, 2, kz).UP.y
    grid:at(kx, 1, kz).UPy, grid:at(kx, 0, kz).UP.y = UPy, UPy
    
    UPz = grid:at(kx, 2, kz).UP.z
    grid:at(kx, 1, kz).UP.z, grid:at(kx, 0, kz).UP.z = UPz, UPz
  end
  end
  
  --X
  for ky = 0, to_y do
  for kz = 0, to_z do
    UPx = grid:at((to_x - 2), ky, kz).UP.x
    grid:at((to_x - 1), ky, kz).UP.x, grid:at(to_x, ky, kz).UP.x = UPx, UPx
    
    UPy = grid:at((to_x - 2), ky, kz).UP.y
    grid:at((to_x - 1), ky, kz).UP.y, grid:at(to_x, ky, kz).UP.y = UPy, UPy
    
    UPz = grid:at((to_x - 2), ky, kz).UP.z
    grid:at((to_x - 1), ky, kz).UP.z, grid:at(to_x, ky, kz).UP.z = UPz, UPz

    UPx = grid:at(2, ky, kz).UP.x
    grid:at(1, ky, kz).UP.x, grid:at(0, ky, kz).UP.x = UPx, UPx
    
    UPy = grid:at(2, ky, kz).UP.y
    grid:at(1, ky, kz).UP.y, grid:at(0, ky, kz).UP.y = UPy, UPy
    
    UPz = grid:at(2, ky, kz).UP.z
    grid:at(1, ky, kz).UP.z, grid:at(0, ky, kz).UP.z = UPz, UPz
  end
  end
end

function P.nonperturbed_UE(grid)
  local to_x = grid:size_x() - 1
  local to_y = grid:size_y() - 1
  local to_z = grid:size_z() - 1
  
  local UEx = 0.0
  local UEy = 0.0
  local UEz = 0.0
  
  -- Z 
  for kx = 0, to_x do
  for ky = 0, to_y do
    UEx = grid:at(kx, ky, (to_z - 3)).UE.x
    grid:at(kx ,ky, (to_z - 2)).UE.x, grid:at(kx ,ky, (to_z - 1)).UE.x, grid:at(kx, ky, to_z).UEx = UEx, UEx, UEx
    
    UEy = grid:at(kx, ky, (to_z - 3)).UE.y
    grid:at(kx ,ky, (to_z - 2)).UE.y, grid:at(kx ,ky, (to_z - 1)).UE.y, grid:at(kx, ky, to_z).UEy = UEy, UEy, UEy
    
    UEz = grid:at(kx, ky, (to_z - 3)).UE.z
    grid:at(kx ,ky, (to_z - 2)).UE.z, grid:at(kx ,ky, (to_z - 1)).UE.z, grid:at(kx, ky, to_z).UE.z = UEz, UEz, UEz
                       
    UEx = grid:at(kx, ky, 2).UE.x
    grid:at(kx, ky, 1).UE.x, grid:at(kx, ky, 0).UE.x = UEx, UEx
    
    UEy = grid:at(kx, ky, 2).UE.y
    grid:at(kx, ky, 1).UE.y, grid:at(kx, ky, 0).UE.y = UEy, UEy
    
    UEz = grid:at(kx, ky, 2).UE.z
    grid:at(kx, ky, 1).UE.z, grid:at(kx, ky, 0).UE.z = UEz, UEz
  end
  end
  
  -- Y
  for kx = 0, to_x do
  for kz = 0, to_z do
    UEx = grid:at(kx, (to_y - 3), kz).UE.x
    grid:at(kx, (to_y - 2), kz).UE.x, grid:at(kx, (to_y - 1), kz).UE.x, grid:at(kx, to_y ,kz).UE.x = UEx, UEx, UEx
    
    UEy = grid:at(kx, (to_y - 3), kz).UE.y
    grid:at(kx, (to_y - 2), kz).UE.y, grid:at(kx, (to_y - 1), kz).UE.y, grid:at(kx, to_y ,kz).UE.y = UEy, UEy, UEy
    
    UEz = grid:at(kx, (to_y - 3), kz).UE.z
    grid:at(kx, (to_y - 2), kz).UE.z, grid:at(kx, (to_y - 1), kz).UE.z, grid:at(kx, to_y ,kz).UE.z = UEz, UEz, UEz
                   
    UEx = grid:at(kx, 2, kz).UE.x
    grid:at(kx, 1, kz).UE.x, grid:at(kx, 0, kz).UE.x = UEx, UEx
    
    UEy = grid:at(kx, 2, kz).UE.y
    grid:at(kx, 1, kz).UE.y, grid:at(kx, 0, kz).UE.y = UEy, UEy
    
    UEz = grid:at(kx, 2, kz).UE.z
    grid:at(kx, 1, kz).UEz, grid:at(kx, 0, kz).UE.z = UEz, UEz
  end
  end
  
  --X
  for ky = 0, to_y do
  for kz = 0, to_z do
    UEx = grid:at(to_x - 2, ky, kz).UE.x
    grid:at((to_x - 1), ky, kz).UE.x, grid:at(to_x, ky, kz).UE.x = UEx, UEx
    
    UEy = grid:at((to_x - 2), ky, kz).UE.y
    grid:at((to_x - 1), ky, kz).UE.y, grid:at(to_x, ky, kz).UE.y = UEy, UEy
    
    UEz = grid:at((to_x - 2), ky, kz).UE.z
    grid:at((to_x - 1), ky, kz).UE.z, grid:at(to_x, ky, kz).UE.z = UEz, UEz

    UEx = grid:at(2, ky, kz).UE.x
    grid:at(1, ky, kz).UE.x, grid:at(0, ky, kz).UE.x = UEx, UEx
    
    UEy = grid:at(2, ky, kz).UE.y
    grid:at(1, ky, kz).UE.y, grid:at(0, ky, kz).UE.y = UEy, UEy
    
    UEz = grid:at(2, ky, kz).UE.z
    grid:at(1, ky, kz).UE.z, grid:at(0, ky, kz).UE.z = UEz, UEz
  end
  end
end

function P.nonperturbed_MF(grid)
  local to_x = grid:size_x() - 1  
  local to_y = grid:size_y() - 1  
  local to_z = grid:size_z() - 1  
                                
  local Bx = 0.0                 
  local By = 0.0                 
  local Bz = 0.0                 

  -- Z 
  for kx = 0, to_x do
  for ky = 0, to_y do
  	Bx = grid:at(kx, ky, (to_z - 3)).B.x
    grid:at(kx ,ky, (to_z - 2)).B.x, grid:at(kx ,ky, (to_z - 1)).B.x, grid:at(kx, ky, to_z).B.x = Bx, Bx, Bx
    
    By = grid:at(kx, ky, (to_z - 3)).B.y
    grid:at(kx ,ky, (to_z - 2)).B.y, grid:at(kx ,ky, (to_z - 1)).B.y, grid:at(kx, ky, to_z).B.y = By, By, By
    
    Bz = grid:at(kx, ky, (to_z - 3)).B.z
    grid:at(kx ,ky, (to_z - 2)).B.z, grid:at(kx ,ky, (to_z - 1)).B.z, grid:at(kx, ky, to_z).B.z = Bz, Bz, Bz
                       
    Bx = grid:at(kx, ky, 2).B.x
    grid:at(kx, ky, 1).B.x, grid:at(kx, ky, 0).B.x = Bx, Bx
    
    By = grid:at(kx, ky, 2).B.y
    grid:at(kx, ky, 1).B.y, grid:at(kx, ky, 0).B.y = By, By
    
    Bz = grid:at(kx, ky, 2).B.z
    grid:at(kx, ky, 1).B.z, grid:at(kx, ky, 0).B.z = Bz, Bz
  end
  end
  
  -- Y
  for kx = 0, to_x do
  for kz = 0, to_z do
    Bx = grid:at(kx, (to_y - 3), kz).B.x
    grid:at(kx, (to_y - 2), kz).B.x, grid:at(kx, (to_y - 1), kz).B.x, grid:at(kx, to_y ,kz).B.x = Bx, Bx, Bx
    
    By = grid:at(kx, (to_y - 3), kz).B.y
    grid:at(kx, (to_y - 2), kz).B.y, grid:at(kx, (to_y - 1), kz).B.y, grid:at(kx, to_y ,kz).B.y = By, By, By
    
    Bz = grid:at(kx, (to_y - 3), kz).B.z
    grid:at(kx, (to_y - 2), kz).B.z, grid:at(kx, (to_y - 1), kz).B.z, grid:at(kx, to_y ,kz).B.z = Bz, Bz, Bz
                   
    Bx = grid:at(kx, 2, kz).B.x
    grid:at(kx, 1, kz).B.x, grid:at(kx, 0, kz).B.x = Bx, Bx
    
    By = grid:at(kx, 2, kz).B.y
    grid:at(kx, 1, kz).B.y, grid:at(kx, 0, kz).B.y = By, By
    
    Bz = grid:at(kx, 2, kz).B.z
    grid:at(kx, 1, kz).B.z, grid:at(kx, 0, kz).B.z = Bz, Bz
  end
  end
  
  --X
  for ky = 0, to_y do
  for kz = 0, to_z do
    Bx = grid:at(to_x - 2, ky, kz).B.x
    grid:at((to_x - 1), ky, kz).B.x, grid:at(to_x, ky, kz).B.x = Bx, Bx
    
    By = grid:at((to_x - 2), ky, kz).B.y
    grid:at((to_x - 1), ky, kz).B.y, grid:at(to_x, ky, kz).B.y = By, By
    
    Bz = grid:at((to_x - 2), ky, kz).B.z
    grid:at((to_x - 1), ky, kz).B.z, grid:at(to_x, ky, kz).B.z = Bz, Bz

    Bx = grid:at(2, ky, kz).B.x
    grid:at(1, ky, kz).B.x, grid:at(0, ky, kz).B.x = Bx, Bx
    
    By = grid:at(2, ky, kz).B.y
    grid:at(1, ky, kz).B.y, grid:at(0, ky, kz).B.y = By, By
    
    Bz = grid:at(2, ky, kz).B.z
    grid:at(1, ky, kz).B.z, grid:at(0, ky, kz).B.z = Bz, Bz
  end
  end
end

function P.nonperturbed_EF(grid)
  local to_x = grid:size_x() - 1  
  local to_y = grid:size_y() - 1  
  local to_z = grid:size_z() - 1  
                                
  local Ex = 0.0                 
  local Ey = 0.0                 
  local Ez = 0.0           
  
  -- Z 
  for kx = 0, to_x do
  for ky = 0, to_y do
  	Ex = grid:at(kx, ky, (to_z - 3)).E.x
  	grid:at(kx, ky, (to_z - 2)).E.x, grid:at(kx, ky, (to_z - 1)).E.x, grid:at(kx, ky, to_z).E.x = Ex, Ex, Ex
    
    Ey = grid:at(kx, ky, (to_z - 3)).E.y
    grid:at(kx, ky, (to_z - 2)).E.y, grid:at(kx, ky, (to_z - 1)).E.y, grid:at(kx, ky, to_z).E.y = Ey, Ey, Ey
    
    Ez = grid:at(kx, ky, (to_z - 3)).E.z
    grid:at(kx, ky, (to_z - 2)).E.z, grid:at(kx, ky, (to_z - 1)).E.z, grid:at(kx, ky, to_z).E.z = Ez, Ez, Ez

    Ex = grid:at(kx, ky, 2).E.x
    grid:at(kx, ky, 1).E.x, grid:at(kx, ky, 0).E.x = Ex, Ex
    
    Ey = grid:at(kx, ky, 2).E.y
    grid:at(kx, ky, 1).E.y, grid:at(kx, ky, 0).E.y = Ey, Ey
    
    Ez = grid:at(kx, ky, 2).E.z
    grid:at(kx, ky, 1).E.z, grid:at(kx, ky, 0).E.z = Ez, Ez
  end
  end
  
  -- Y
  for kx = 0, to_x do
  for kz = 0, to_z do
    Ex = grid:at(kx, (to_y - 3), kz).E.x
    grid:at(kx, (to_y - 2), kz).E.x, grid:at(kx, (to_y - 1), kz).E.x, grid:at(kx, to_y, kz).E.x = Ex, Ex, Ex
    
    Ey = grid:at(kx, (to_y - 3), kz).E.y
    grid:at(kx, (to_y - 2), kz).E.y, grid:at(kx, (to_y - 1), kz).E.y, grid:at(kx, to_y, kz).E.y = Ey, Ey, Ey
    
    Ez = grid:at(kx, (to_y - 3), kz).E.z
    grid:at(kx, (to_y - 2), kz).E.z, grid:at(kx, (to_y - 1), kz).E.z, grid:at(kx, to_y, kz).E.z = Ez, Ez, Ez
                   
    Ex = grid:at(kx, 2, kz).E.x
    grid:at(kx, 1, kz).E.x, grid:at(kx, 0, kz).E.x = Ex, Ex
    
    Ey = grid:at(kx, 2, kz).E.y
    grid:at(kx, 1, kz).E.y, grid:at(kx, 0, kz).E.y = Ey, Ey
    
    Ez = grid:at(kx, 2, kz).E.z
    grid:at(kx, 1, kz).E.z, grid:at(kx, 0, kz).E.z = Ez, Ez
  end
  end
  
  --X
  for ky = 0, to_y do
  for kz = 0, to_z do
    Ex = grid:at((to_x - 2), ky, kz).E.x
    grid:at((to_x - 1), ky, kz).E.x, grid:at(to_x, ky, kz).E.x = Ex, Ex
    
    Ey = grid:at((to_x - 2), ky, kz).E.y
    grid:at((to_x - 1), ky, kz).E.y, grid:at(to_x, ky, kz).E.y = Ey, Ey
    
    Ez = grid:at((to_x - 2), ky, kz).E.z
    grid:at((to_x - 1), ky, kz).E.z, grid:at(to_x, ky, kz).E.z = Ez, Ez

    Ex = grid:at(2, ky, kz).E.x
    grid:at(1, ky, kz).E.x, grid:at(0, ky, kz).E.x = Ex, Ex
    
    Ey = grid:at(2, ky, kz).E.y
    grid:at(1, ky, kz).E.y, grid:at(0, ky, kz).E.y = Ey, Ey
    
    Ez = grid:at(2, ky, kz).E.z
    grid:at(1, ky, kz).E.z, grid:at(0, ky, kz).E.z = Ez, Ez

  end
  end
end
