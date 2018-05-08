--[[
NAME
  print_mpi

FUNCTIONS
  print_root
  fprint_root
    
NOTES
    Package for print to standard output (stdout) and file
]]

local P = {}

if _REQUIREDNAME == nil then
    print_mpi = P
else
    _G[_REQUIREDNAME] = P
end

function P.print_root(proc_idx, ...)
  if proc_idx == 0 then 
    for i,v in ipairs(arg) do
      io.write(v)
    end
  end
end
   
function P.fprint_root(proc_idx, out_file, ...)
  if proc_idx == 0 then 
    for i,v in ipairs(arg) do
       out_file:write(v)
    end
  end
end