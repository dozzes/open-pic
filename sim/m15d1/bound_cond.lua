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
  Non-perturbed (zero-gradient) boundary conditions for the staggered C-grid.
  Every boundary node, including edges and corners, is written exactly once.
]]

local P = {}

if _REQUIREDNAME == nil then
  bound_cond = P
else
  _G[_REQUIREDNAME] = P
end

local function clamp_source(index, upper_source)
  if index < 2 then
    return 2
  elseif index > upper_source then
    return upper_source
  end
  return index
end

-- Copy one scalar component to the boundary. The upper source index depends
-- on the component's C-grid centering. X owns edges and corners, Y excludes
-- the X boundary, and Z excludes both X and Y, preventing overwrites.
local function copy_component(grid, upper_x, upper_y, upper_z, reader, writer)
  local to_x = grid:size_x() - 1
  local to_y = grid:size_y() - 1
  local to_z = grid:size_z() - 1

  local function copy_x(x)
    local sx = clamp_source(x, upper_x)
    for y = 0, to_y do
    for z = 0, to_z do
      local sy = clamp_source(y, upper_y)
      local sz = clamp_source(z, upper_z)
      writer(grid:at(x, y, z), reader(grid:at(sx, sy, sz)))
    end
    end
  end

  for x = 0, 1 do copy_x(x) end
  for x = upper_x + 1, to_x do copy_x(x) end

  local function copy_y(y)
    local sy = clamp_source(y, upper_y)
    for x = 2, upper_x do
    for z = 0, to_z do
      local sz = clamp_source(z, upper_z)
      writer(grid:at(x, y, z), reader(grid:at(x, sy, sz)))
    end
    end
  end

  for y = 0, 1 do copy_y(y) end
  for y = upper_y + 1, to_y do copy_y(y) end

  local function copy_z(z)
    local sz = clamp_source(z, upper_z)
    for x = 2, upper_x do
    for y = 2, upper_y do
      writer(grid:at(x, y, z), reader(grid:at(x, y, sz)))
    end
    end
  end

  for z = 0, 1 do copy_z(z) end
  for z = upper_z + 1, to_z do copy_z(z) end
end

local function copy_face_vector(grid, field)
  local to_x = grid:size_x() - 1
  local to_y = grid:size_y() - 1
  local to_z = grid:size_z() - 1

  -- X component: [i, j + 1/2, k + 1/2]
  copy_component(grid, to_x - 2, to_y - 3, to_z - 3,
    function(cell) return cell[field].x end,
    function(cell, value) cell[field].x = value end)

  -- Y component: [i + 1/2, j, k + 1/2]
  copy_component(grid, to_x - 3, to_y - 2, to_z - 3,
    function(cell) return cell[field].y end,
    function(cell, value) cell[field].y = value end)

  -- Z component: [i + 1/2, j + 1/2, k]
  copy_component(grid, to_x - 3, to_y - 3, to_z - 2,
    function(cell) return cell[field].z end,
    function(cell, value) cell[field].z = value end)
end

function P.nonperturbed_NP(grid)
  local to_x = grid:size_x() - 1
  local to_y = grid:size_y() - 1
  local to_z = grid:size_z() - 1

  -- Cell centered: [i + 1/2, j + 1/2, k + 1/2]
  copy_component(grid, to_x - 3, to_y - 3, to_z - 3,
    function(cell) return cell.NP end,
    function(cell, value) cell.NP = value end)
end

function P.nonperturbed_UP(grid)
  copy_face_vector(grid, "UP")
end

function P.nonperturbed_UE(grid)
  copy_face_vector(grid, "UE")
end

function P.nonperturbed_EF(grid)
  copy_face_vector(grid, "E")
end

function P.nonperturbed_MF(grid)
  local to_x = grid:size_x() - 1
  local to_y = grid:size_y() - 1
  local to_z = grid:size_z() - 1

  -- Bx: [i + 1/2, j, k]
  copy_component(grid, to_x - 3, to_y - 2, to_z - 2,
    function(cell) return cell.B.x end,
    function(cell, value) cell.B.x = value end)

  -- By: [i, j + 1/2, k]
  copy_component(grid, to_x - 2, to_y - 3, to_z - 2,
    function(cell) return cell.B.y end,
    function(cell, value) cell.B.y = value end)

  -- Bz: [i, j, k + 1/2]
  copy_component(grid, to_x - 2, to_y - 2, to_z - 3,
    function(cell) return cell.B.z end,
    function(cell, value) cell.B.z = value end)
end

return P
