--[[
NAME
  parts_count

FUNCTIONS
  get_cloud_parts_num(R0, dist)
  get_backgr_parts_num(grid, Bz0, ratio_Bmax_to_Bmin, center_x, center_y, center_z)
  get_cyl_parts_num(R0, H, dist)

NOTES
    Package for particles initialization
]]

local P = {}

if _REQUIREDNAME == nil then
    parts_count = P
else
    _G[_REQUIREDNAME] = P
end

function P.get_cloud_parts_num(R0, dist)
  local to_R = R0 + dist
  local parts_num = 0
  for x = 0.0, to_R, dist  do
  for y = 0.0, to_R, dist  do
  for z = 0.0, to_R, dist  do
    local r = math.sqrt(x*x + y*y + z*z)
    if (r <= cloud_R0) then 
      parts_num = parts_num + 1 -- // +++ octant
      if z ~= 0.0 then parts_num = parts_num + 1 end -- // ++- octant
      if y ~= 0.0 then parts_num = parts_num + 1 end -- // +-+ octant
      if x ~= 0.0 then parts_num = parts_num + 1 end -- // -++ octant
      if x ~= 0.0 and y ~= 0 then parts_num = parts_num + 1 end -- // --+ octant
      if x ~= 0.0 and z ~= 0 then parts_num = parts_num + 1 end -- // -+- octant
      if y ~= 0.0 and z ~= 0 then parts_num = parts_num + 1 end -- // +-- octant
      if x ~= 0 and y ~= 0.0 and z ~= 0 then parts_num = parts_num + 1 end -- // --- octant
    end
  end
  end
  end
  return parts_num
end  -- cloud_parts_num

 function P.get_backgr_parts_num(h, length_x, length_y, length_z, center_x, center_y, center_z, cloud_R0, dist)
   local to_x = 0.5*length_x - h
   local to_y = 0.5*length_y - h
   local to_z = 0.5*length_z - h
  
   local parts_num = 0
   for x = 0, to_x, dist do
   for y = 0, to_y, dist do
   for z = 0, to_z, dist do
     cx = center_x - x
     cy = center_y - y
     cz = center_z - z
     if math.sqrt(cx*cx + cy*cy + cz*cz) > cloud_R0 then
       parts_num = parts_num + 1 -- // +++ octant
       if z ~= 0.0 then parts_num = parts_num + 1 end -- // ++- octant
       if y ~= 0.0 then parts_num = parts_num + 1 end -- // +-+ octant
       if x ~= 0.0 then parts_num = parts_num + 1 end -- // -++ octant
       if x ~= 0.0 and y ~= 0 then parts_num = parts_num + 1 end -- // --+ octant
       if x ~= 0.0 and z ~= 0 then parts_num = parts_num + 1 end -- // -+- octant
       if y ~= 0.0 and z ~= 0 then parts_num = parts_num + 1 end -- // +-- octant
       if x ~= 0 and y ~= 0.0 and z ~= 0 then parts_num = parts_num + 1 end -- // --- octant
     end
   end  
   end
   end
   return parts_num
end  -- backgr_parts_num

