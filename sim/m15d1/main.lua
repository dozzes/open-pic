  --[[

    Main configuration script for opic3d
    Defines task for plasma cloud expansion in magnetized backround.

    All opic3d variable names have f prefix "pic_":
      pic_parameters - common parameters for opic3d;
      pic_grid - computational grid;
      pic_particles - particles array;
      pic_particle_groups - particle groups;
      pic_marker_particles - marked particles. Used for particle trajectory analysis.

   The detailed description is availble at  http://comphys.narod.ru/open_pic_3d/open_pic_3d.htm.

  ]]

  require("print_mpi")

  function on_iteration_begin() -- called by opic3d.exe in the beginning of each iteration on time
  end

  function on_iteration_end() -- called in the end of each iteration on time
    local t = pic_parameters.current_time_step
    local st = pic_parameters.save_time_steps

    if (t == 1) then
      os.execute("mkdir diag")
    end

    local ts = (t % st)

    if (ts == 0) or (t == st) then
      os.execute("cp *markers*.dat ./diag")
      os.execute("cp at*.dat ./diag")
      os.execute("cp *_x_" .. cloud_x_node .. "_grd_*.dat ./diag")
      os.execute("cp *_y_" .. cloud_y_node .. "_grd_*.dat ./diag")
      os.execute("cp *_z_" .. cloud_z_node .. "_grd_*.dat ./diag")
      io.write("backup diagnostic data...\n")
      os.execute("tar --remove-files -czvf  " .. t .. ".tar.gz *.dat")
    end


    if (t == 1) or (ts == 0) or (t == st) then
      local wmf = 0.0
      local wel = 0.0
      for kx = 0, (pic_grid:size_x() - 1) do
      for ky = 0, (pic_grid:size_y() - 1) do
      for kz = 0, (pic_grid:size_z() - 1) do
      	local b = pic_grid:at(kx,ky,kz).B:abs()
        wmf = wmf + b*b

        local ue = pic_grid:at(kx,ky,kz).UE:abs()
        wel = wel + ue*ue*pic_grid:at(kx,ky,kz).NP
      end
      end
      end
      wmf = 0.125*wmf/pi/cloud_W0
      wel = 0.5*me*wel/cloud_W0
      print_mpi.print_root(proc_idx, "\nwe = ", wel, "\n")

      local wc = 0.0
      local wb = 0.0
      for k = 0, (pic_particles.size-1) do
      	p = pic_particles:at(k)
      	v_abs = p.v:abs()
      	if (p.grp == "cloud") then
      		wc = wc + v_abs*v_abs
        elseif (p.grp == "backgr") then
          wb = wb + v_abs*v_abs
        end
      end
      wc = 0.5*wc*cloud_mi*ni/cloud_W0
      wb = 0.5*wb*backgr_mi*ni/cloud_W0

      local f = io.open("energy.txt", "a")

      if (t == 1) then
      	 f:write("t\tW_cloud\tW_backgr\tW_mf\tW_e\n")
      end

      local tnd = t * pic_parameters.tau/pic_parameters.T_scale

      f:write(tnd, "\t", wc, "\t", wb, "\t", wmf, "\t", wel, "\n")
      f:flush()
    end
  end

  -- pos - spatial coordinate as DbLVector
  function at_45(pos) -- called grid save filter at x=y=z by opic3d.exe to set boundary coditions for density (NP)
    return (pos.x == pos.y and pos.y == pos.z)
  end

  function on_set_boundary_NP() -- called by opic3d.exe to set boundary coditions for density (NP)
    require("bound_cond")
    bound_cond.nonperturbed_NP(pic_grid)
  end

  function on_set_boundary_UP() -- called by opic3d.exe to set boundary coditions for ion velocity (UP)
  	require("bound_cond")
    bound_cond.nonperturbed_UP(pic_grid)
  end

  function on_set_boundary_UE() -- called by opic3d.exe to set boundary coditions for electron velocity (UE)
  	require("bound_cond")
    bound_cond.nonperturbed_UE(pic_grid)
  end

  function on_set_boundary_EF() -- called by opic3d.exe to set boundary coditions for electic field (EF)
  	require("bound_cond")
    bound_cond.nonperturbed_EF(pic_grid)
  end

  function on_set_boundary_MF() -- called by opic3d.exe to set boundary coditions for magnetic field (MF)
     require("bound_cond")
     bound_cond.nonperturbed_MF(pic_grid)
  end

  function rnd(part_dist)
     local w = 0.5*part_dist*(0.5 - math.random())
     return w
  end

  function rnd_r(r, part_dist)
     local w = r + 0.5*part_dist*(0.5 - math.random())
    return w
  end

do
  proc_idx = pic_parameters.process_idx;
  print_mpi.print_root(proc_idx, "\n*** start ***\n\n")

  -- physical constants
  c  = 2.9979e+10 -- light velocity for vacuum (cm/sec)
  e  = 4.8032e-10   -- electron charge
  me = 9.1e-28     -- electron charge
  mp = 1.6726e-24 -- proton mass

  pi = math.pi

  -- background parameters
  backgr_mi   = 1 * mp  -- background ions mass (in proton mass units)
  backgr_dens = 1.0e+14 -- background ions density (cm^-3)
  backgr_Z    = 1       -- background ions	charge (in proton charge units)

  -- magnetic field
  Bz0 = 100.0 -- initial magnetic field (Gs)
  B0  = Bz0   -- initial maximum of magnetic field
  V_A0 = B0 / math.sqrt(4*pi*backgr_dens*backgr_mi) -- Alfven velocity in background

  -- cloud parameters
  cloud_mi       = 1 * mp    -- cloud ions	mass (in proton	mass units)
  cloud_ions_num = 1.66E+19  -- cloud	ions number
  cloud_Z        = 1         -- cloud ions	charge (in proton charge units)
  cloud_v_min    = 0.0       -- cloud min velocity (cm/sec)
  Ma             = 15.0      -- Mach-Alfven number
  cloud_v_max    = Ma * V_A0 -- cloud max velocity (cm/sec)

  -- simulation parameters
  V_max = math.max(cloud_v_max, V_A0) -- max velocity
  break_times          = 4.9  -- number of cloud braking times to simulate
  points_on_RL         = 10   -- number of points on Larmor radiis
  h_R0                 = 0.1  -- initial cloud radius (in units of h)
  cloud_parts_on_step  = 1000 -- simulation particle for cloud per spatial step
  backgr_parts_on_step = 2    -- simulation particle for background per spatial step

  -- spatial step (h = hx = hy = hz)
  RL = cloud_mi * cloud_v_max * c / (cloud_Z * e * B0) -- Larmour radius for cloud ions
  backgr_Wpi   = backgr_Z * e * 2 * math.sqrt(pi * backgr_dens / backgr_mi) -- background ion plasma frequency
  backgr_c_Wpi = c / backgr_Wpi -- background ion inertial length
  backgr_Wci   = backgr_Z * e * B0 / (backgr_mi * c) -- background ion cyclotron frequency
  RL_b = cloud_v_max/backgr_Wci -- Larmour radius for background ions
  h = math.max(0.0*backgr_c_Wpi, RL / points_on_RL) -- spatial step of grid (cm)
  pic_grid.step = h
  print_mpi.print_root(proc_idx, "h = ", h, "\n")
  print_mpi.print_root(proc_idx, "backgr_c_Wpi = ", backgr_c_Wpi, "\n")

  cloud_R0 = h * h_R0 -- initial cloud radius (cm)
  cell_volume  = h*h*h -- cell volume
  cloud_volume = 4.0 / 3.0 * pi * math.pow(cloud_R0, 3) -- cloud volume (cm^3)
  cloud_dens   = cloud_ions_num / cloud_volume -- cloud dnsity (cm^(-3))
  cloud_Wpi    = 2 * cloud_Z * e * math.sqrt(pi * cloud_dens / cloud_mi) -- cloud ion plasma frequency
  cloud_c_Wpi  = c / cloud_Wpi -- cloud ion inertial length
  cloud_Wci    = cloud_Z * e * B0/(cloud_mi * c) -- cloud ion cyclotron frequency

  -- timestep restricted by the local CFL condition on plasma wave propagation
  Wci_max = math.max(cloud_Wci, backgr_Wci) -- max ion cyclotron frequency
  CFL_tau = 0.5 * h / V_max -- timestep from linear CLF condition
  Wci_tau = 0.1 / Wci_max   -- timestep must be < Wci_max^(-1)
  CFL_tau_loc = math.pow(h / math.max(cloud_c_Wpi, backgr_c_Wpi), 2) / (pi * math.sqrt(3) * Wci_max)
  tau = 0.2 * math.min(CFL_tau, Wci_tau, CFL_tau_loc) -- best tau

  print_mpi.print_root(proc_idx, "CFL_tau = ", CFL_tau, "\n")
  print_mpi.print_root(proc_idx, "Wci_tau = ", Wci_tau, "\n")
  print_mpi.print_root(proc_idx, "CFL_tau_loc = ", CFL_tau_loc, "\n")
  print_mpi.print_root(proc_idx, "tau = ", tau, "\n")

  pic_parameters.tau = tau -- set pic simulation timestemp

  cloud_part_dist = h / cloud_parts_on_step -- distance between nearest cloud "macroparticle"

  require("parts_count") -- require cloud_parts package
  cloud_parts_num = parts_count.get_cloud_parts_num(cloud_R0, cloud_part_dist) -- necessary "macroparticles" count

  cloud_mass = cloud_ions_num * cloud_mi -- cloud mass
  cloud_W0 = 0.3 * cloud_mass * V_max * V_max -- cloud energy

  Rb    = math.pow(6 * cloud_W0 / (B0 * B0), 1/3) -- braking on magnetic field radius
  Rb1   = math.pow(3 * cloud_ions_num / (4 * pi * backgr_dens * backgr_Z), 1.0/3.0) -- braking on magnetic field radius
  Rg    = Rb / math.pow(Ma, 2.0/3.0) -- gasdynamic braking radius
  Rt    = math.min(Rb, Rg) -- actual braking radius
  Tb    = Rb / V_max -- braking time on magnetic field
  Tg    = Rg / V_max -- gasdynamic braking time
  Tt    = math.min(Tb, Tg)  -- actual braking time
  Delta = math.pow(Rg / RL, 2) -- magtetolaminar parameter
  Eps   = RL / Rb -- epsilon parameter

  --- print basic parameters
  print_mpi.print_root(proc_idx, "RL = ", RL, "\n")
  print_mpi.print_root(proc_idx, "Rb = ", Rb, "\n")
  print_mpi.print_root(proc_idx, "Rb1 = ", Rb1, "\n")
  print_mpi.print_root(proc_idx, "Delta = ", Delta, "\n")
  print_mpi.print_root(proc_idx, "backgr_c_Wpi = ", backgr_c_Wpi, "\n")
  print_mpi.print_root(proc_idx, "cloud_c_Wpi = ", cloud_c_Wpi, "\n")

  -- set PIC simulation timestemps number
  sim_times = 2.0
  pic_parameters.time_steps = sim_times * break_times * Tt / tau
  print_mpi.print_root(proc_idx, "time_steps = ", pic_parameters.time_steps, "\n")

  -- set PIC simulation data save timestemps period
  pic_parameters.save_time_steps = 10

  -- Wci*tau < 0.2
  -- h > c/Wpi
  -- save parameters
  local f = io.open("lua_params.txt", "w")
  f:write("h                = ", h, "(shall be  > ", cloud_c_Wpi, ", backgr_c_Wpi = ", backgr_c_Wpi, ")\n")
  f:write("cloud_c_Wpi      = ", cloud_c_Wpi, "\n")
  f:write("cloud_Wpi        = ", cloud_Wpi, "\n")
  f:write("cloud_Wci        = ", cloud_Wci, "\n")
  f:write("cloud_c_Wci      = ", c/cloud_Wci, "\n")
  f:write("cloud_Wci*tau    = ", (cloud_Wci*tau), " ( shall be < 0.2 )\n")
  f:write("cloud_c_Wci      = ", (c/cloud_Wci), "\n")
  f:write("cloud_R0         = ", cloud_R0, "\n")
  f:write("cloud_W0         = ", cloud_W0, "\n")
  f:write("cloud_part_dist  = ", cloud_part_dist, "\n")
  f:write("cloud_parts_num  = ", cloud_parts_num, "\n")
  local backgr_Wpe = e * 2 * math.sqrt(pi * backgr_dens / me) -- background electron plasma frequency
  f:write("backgr_Wpe       = ", backgr_Wpe, "\n")
  f:write("backgr_c_Wpe     = ", c/backgr_Wpe, "\n")
  f:write("backgr_Wpi       = ", backgr_Wpi, "\n")
  f:write("backgr_c_Wpi     = ", backgr_c_Wpi, "\n")
  f:write("backgr_Wci       = ", backgr_Wci, "\n")
  f:write("backgr_Wci*tau   = ", (backgr_Wci*tau), "\n")
  f:write("RL               = ", RL, "\n")
  f:write("RL_b             = ", RL_b, "\n")
  f:write("h                = ", h, "\n")
  f:write("tau              = ", tau, "\n")
  f:write("CFL_tau          = ", CFL_tau, "\n")
  f:write("Wci_tau          = ", Wci_tau, "\n")
  f:write("CFL_tau_loc      = ", CFL_tau_loc, "\n")
  f:write("Rb               = ", Rb, "\n")
  f:write("Rb1              = ", Rb1, "\n")
  f:write("Rg               = ", Rg, "\n")
  f:write("Rt               = ", Rt, "\n")
  f:write("Tb               = ", Tb, "\n")
  f:write("Tg               = ", Tg, "\n")
  f:write("Tt               = ", Tt, "\n")
  f:write("Ma               = ", Ma, "\n")
  f:write("V_A0             = ", V_A0, "\n")
  f:write("Delta            = ", Delta, "\n")
  f:write("Eps              = ", Eps, "\n")
  f:flush()

  -- calculate simulation grid size
  x_nodes = math.floor(2 * V_max * break_times * Tt / h + 2)
  if x_nodes % 2 == 0 then x_nodes = x_nodes + 1 end
  y_nodes = x_nodes
  z_nodes = x_nodes

  f:write("x_nodes = ", x_nodes, "\n")
  f:write("y_nodes = ", y_nodes, "\n")
  f:write("z_nodes = ", z_nodes, "\n")

  print_mpi.print_root(proc_idx, "x_nodes = ", x_nodes, "\n")
  print_mpi.print_root(proc_idx, "y_nodes = ", y_nodes, "\n")
  print_mpi.print_root(proc_idx, "z_nodes = ", z_nodes, "\n")

  -- set grid size
  pic_grid:resize(x_nodes, y_nodes, z_nodes)
  pic_grid:set_boundary_state(Cell.cs_absorptive)

  -- calculate cloud center position
  length_x = pic_grid:length_x()
  length_y = pic_grid:length_y()
  length_z = pic_grid:length_z()

  cloud_x_node = math.floor(0.5 * x_nodes)
  cloud_y_node = math.floor(0.5 * y_nodes)
  cloud_z_node = math.floor(0.5 * z_nodes)

  cloud_x = h * (cloud_x_node)
  cloud_y = h * (cloud_y_node)
  cloud_z = h * (cloud_z_node)

  f:write("length_x = ", length_x, "\n")
  f:write("length_y = ", length_y, "\n")
  f:write("length_z = ", length_z, "\n")
  f:write("cloud_x = ", cloud_x, "\n")
  f:write("cloud_y = ", cloud_y, "\n")
  f:write("cloud_z = ", cloud_z, "\n")
  f:write("cloud_x_node = ", cloud_x_node, "\n")
  f:write("cloud_y_node = ", cloud_y_node, "\n")
  f:write("cloud_z_node = ", cloud_z_node, "\n")

  backgr_part_dist = h / backgr_parts_on_step -- distance between nearest background "macroparticle"

  -- background "macroparticles" count
  backgr_parts_num = parts_count.get_backgr_parts_num(h, length_x, length_y, length_z,
                                                      cloud_x, cloud_y, cloud_z, cloud_R0,
                                                      backgr_part_dist)

  f:write("backgr_parts_num = ", backgr_parts_num, "\n")

  -- total "macroparticles" count
  total_parts_num = cloud_parts_num + backgr_parts_num
  print_mpi.print_root(proc_idx, "cloud_parts_num = ", cloud_parts_num, "\n")
  print_mpi.print_root(proc_idx, "backgr_parts_num = ", backgr_parts_num, "\n")
  print_mpi.print_root(proc_idx, "total_parts_num = ", total_parts_num, "\n")

  -- check if memory available (change for you purpose)
  mem_bytes_needed = (total_parts_num * 96 + x_nodes * y_nodes * z_nodes * 80)
  print_mpi.print_root(proc_idx, (math.ceil(mem_bytes_needed/1024/1024)), " MB needed.\n")
  if (mem_bytes_needed > 2e+9)then
    print_mpi.print_root(proc_idx, "\nMemory allocation fail is possible. Proceed anyway? [y/n] ")
    yn = io.read()
    if (yn ~= "y") and (yn ~= "Y") then return end
  end

  f:write("total_parts_num  = ", total_parts_num, "\n")

  -- set pic simulation "macroparticles" count
  pic_particles.size = total_parts_num

  -- initialize cloud particles
  io.write("\n*** initialize cloud particles ***\n")

  ni = cloud_ions_num / cloud_parts_num
  f:write("cloud_ni = ", ni, "\n")

  pic_particle_groups:create_group("cloud", 1.0, 1.0, ParticleGroups.save_grid_values)
  to_R = cloud_R0 + cloud_part_dist
  cloud_W0_num = 0
  pk = 0

  for x = 0.0, to_R, cloud_part_dist do
  for y = 0.0, to_R, cloud_part_dist do
  for z = 0.0, to_R, cloud_part_dist do

    local r = math.sqrt(x*x + y*y + z*z)

    if (r <= cloud_R0) then

    	local px = x + rnd(cloud_part_dist)
      local py = y + rnd(cloud_part_dist)
      local pz = z + rnd(cloud_part_dist)

      vx = cloud_v_min + cloud_v_max * px / cloud_R0
      vy = cloud_v_min + cloud_v_max * py / cloud_R0
      vz = cloud_v_min + cloud_v_max * pz / cloud_R0

      -- octant 000
      r = DblVector((cloud_x + px), (cloud_y + py), (cloud_z + pz))
      v = DblVector(vx, vy, vz)
      pic_particles:set(pk, "cloud", r, v, ni)
      if (z == 0.0) then
        pic_marker_particles:add("cloud_z=0", pk)
        if (x == 0.0) or (y == 0.0) or (x == y) then pic_marker_particles:add("cloud_mrk", pk) end
      end
      pk = pk + 1

      -- octant 001
      if (z ~= 0.0) then
      	r = DblVector((cloud_x + px), (cloud_y + py), (cloud_z - pz))
      	v = DblVector(vx, vy, (-vz))
        pic_particles:set(pk, "cloud", r, v, ni)
        pk = pk + 1
      end

      -- octant 010
      if (y ~= 0.0) then
      	r = DblVector((cloud_x + px), (cloud_y - py), (cloud_z + pz))
      	v = DblVector(vx, (-vy), vz)
        pic_particles:set(pk, "cloud", r, v, ni)
        if (z == 0.0) then
          pic_marker_particles:add("cloud_z=0", pk)
          if (x == 0.0) or (x == y) then pic_marker_particles:add("cloud_mrk", pk) end
        end
        pk = pk + 1
      end

      -- octant 011
      if (y ~= 0.0) and (z ~= 0.0) then
      	r = DblVector((cloud_x + px), (cloud_y - py), (cloud_z - pz))
      	v = DblVector(vx, (-vy), (-vz))
        pic_particles:set(pk, "cloud", r, v, ni)
        pk = pk + 1
      end

      -- octant 100
      if (x ~= 0.0) then
      	r = DblVector((cloud_x - px), (cloud_y + py), (cloud_z + pz))
      	v = DblVector((-vx), vy, vz)
        pic_particles:set(pk, "cloud", r, v, ni)
        if (z == 0.0) then
           pic_marker_particles:add("cloud_z=0", pk)
           if (y == 0.0) or (x == y) then pic_marker_particles:add("cloud_mrk", pk) end
        end
        pk = pk + 1
      end

      -- octant 101
      if (x ~= 0.0) and (z ~= 0.0) then
      	r = DblVector((cloud_x - px), (cloud_y + py), (cloud_z - pz))
      	v = DblVector((-vx), (vy), (-vz))
        pic_particles:set(pk, "cloud", r, v, ni)
        pk = pk + 1
      end

      -- octant 110
      if (x ~= 0.0) and (y ~= 0.0) then
      	r = DblVector((cloud_x - px), (cloud_y - py), (cloud_z + pz))
      	v = DblVector((-vx), (-vy), vz)
        pic_particles:set(pk, "cloud", r, v, ni)
        if (z == 0.0) then
          pic_marker_particles:add("cloud_z=0", pk)
          if (x == y) then pic_marker_particles:add("cloud_mrk", pk) end
        end
        pk = pk + 1
      end

      -- octant 111
      if (x ~= 0.0) and (y ~= 0.0) and (z ~= 0.0) then
      	r = DblVector((cloud_x - px), (cloud_y - py), (cloud_z - pz))
      	v = DblVector((-vx), (-vy), (-vz))
        pic_particles:set(pk, "cloud", r, v, ni)
        pk = pk + 1
      end

      v_abs = v:abs();
      cloud_W0_num = cloud_W0_num + v_abs*v_abs;

    end -- if (r <= cloud_R0)
  end
  end
  end

  print_mpi.print_root(proc_idx, "pk = ", pk, "\n")

  cloud_W0_num = 4*cloud_W0_num*ni*cloud_mi
  print_mpi.print_root(proc_idx, "cloud_W0_num = ", cloud_W0_num, "\n")
  print_mpi.print_root(proc_idx, "cloud_W0 = ", cloud_W0, "\n")

  f:write("cloud_W0 = ", cloud_W0, "\n")
  f:write("cloud_W0_num   = ", cloud_W0_num, "\n")

  print_mpi.print_root(proc_idx, "\n*** initialize background particles ***\n")

  -- initialize background particles
  backgr_ni = backgr_dens * cell_volume / math.pow(backgr_parts_on_step, 3)
  pic_parameters.dens_cutoff = 0.00001* backgr_ni
  f:write("backgr_ni = ", backgr_ni, "\n")
  print_mpi.print_root(proc_idx, "dens_cutoff = ", pic_parameters.dens_cutoff, "\n")

  pic_particle_groups:create_group("backgr", 1.0, 1.0, ParticleGroups.save_grid_values)

  backgr_vx = 0.0
  backgr_v = DblVector(backgr_vx, 0.0, 0.0)
  pk = cloud_parts_num
  to_x = 0.5*length_x - h
  to_y = 0.5*length_y - h
  to_z = 0.5*length_z - h

  print_mpi.print_root(proc_idx, "pk = ", pk, "\n")

  for x = 0, to_x, backgr_part_dist do
  for y = 0, to_y, backgr_part_dist do
  for z = 0, to_z, backgr_part_dist do

     local cx = cloud_x - x
     local cy = cloud_y - y
     local cz = cloud_z - z

     if (math.sqrt(cx*cx + cy*cy + cz*cz) > cloud_R0) then

       -- octant 000
       local px = x --+ rnd(backgr_part_dist)
       local py = y --+ rnd(backgr_part_dist)
       local pz = z --+ rnd(backgr_part_dist)

       local r = DblVector((cloud_x + px), (cloud_y + py), (cloud_z + pz))
       pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
       if (z == 0.0) then
         pic_marker_particles:add("backgr_z=0", pk)
         if (x == 0.0) or (y == 0) or (x == y) then pic_marker_particles:add("backgr_mrk", pk) end
       end
       pk = pk + 1

       -- octant 001
       if (z ~= 0.0) then
         r = DblVector((cloud_x + px), (cloud_y + py), (cloud_z - pz))
         pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
         pk = pk + 1
       end

       -- octant 010
       if (y ~= 0.0) then
         r = DblVector((cloud_x + px), (cloud_y - py), (cloud_z + pz))
         pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
         if (z == 0.0) then
           pic_marker_particles:add("backgr_z=0", pk)
           if (x == 0.0) or (x == y) then pic_marker_particles:add("backgr_mrk", pk) end
         end
         pk = pk + 1
       end

       -- octant 011
       if (y ~= 0.0) and (z ~= 0.0) then
         r = DblVector((cloud_x + px), (cloud_y - py), (cloud_z - pz))
         pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
         pk = pk + 1
       end

       -- octant 100
       if (x ~= 0.0) then
         r = DblVector((cloud_x - px), (cloud_y + py), (cloud_z + pz))
         pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
         if (z == 0.0) then
           pic_marker_particles:add("backgr_z=0", pk)
           if (y == 0.0) or (x == y) then pic_marker_particles:add("backgr_mrk", pk) end
         end
         pk = pk + 1
       end

       -- octant 101
       if (x ~= 0.0) and (z ~= 0.0) then
         r = DblVector((cloud_x - px), (cloud_y + py), (cloud_z - pz))
         pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
         pk = pk + 1
       end

       -- octant 110
       if (x ~= 0.0) and (y ~= 0.0) then
         r = DblVector((cloud_x - px), (cloud_y - py), (cloud_z + pz))
         pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
         if (z == 0.0) then
           pic_marker_particles:add("backgr_z=0", pk)
           if (x == y) then pic_marker_particles:add("backgr_mrk", pk) end
         end
         pk = pk + 1
       end

       -- octant 111
       if x ~= 0.0 and y ~= 0.0 and z ~= 0.0 then
         r = DblVector((cloud_x - px), (cloud_y - py), (cloud_z - pz))
         pic_particles:set(pk, "backgr", r, backgr_v, backgr_ni)
         pk = pk + 1
       end
    end
  end
  end
  end
  print_mpi.print_root(proc_idx, "backround particles = ", (pk - cloud_parts_num), "\n")

  -- initialize EM field
  print_mpi.print_root(proc_idx, "\n*** initialize EM field ***\n")

  require("em_field")
  --em_field.uniform(pic_grid, 0.0, 0.0, 0.0, 0.0, 0.0, Bz0)
  em_field.dipole_equator(pic_grid, Bz0, 100.0)


  pic_parameters.save_all_particles = false
  pic_parameters.save_whole_grid = true

  pic_parameters.save_grid_x_plains = true
  pic_parameters.save_grid_y_plains = true
  pic_parameters.save_grid_z_plains = true

  pic_parameters.L_scale = backgr_c_Wpi
  pic_parameters.T_scale = 1.0/backgr_Wci
  pic_parameters.U_scale = V_A0
  pic_parameters.N_scale = backgr_dens
  pic_parameters.E_scale = (V_A0/c)*B0
  pic_parameters.B_scale = B0

  --pic_grid_save_filters:add("at_45")

  print_mpi.print_root(proc_idx, "\n*** end ***\n\n")
  io.flush()
  f:flush()
end
