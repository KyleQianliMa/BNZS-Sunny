using Sunny
using LinearAlgebra
using GLMakie, Plots

############## load cif file ##############
cif_file = joinpath(@__DIR__, "BNZS_SSL.cif")
cryst = Crystal(cif_file, symprec=1e-3)
Nd_subcryst = subcrystal(cryst, "Nd1");
# view_crystal(Nd_subcryst, 7.0)

dim = (1, 1, 1)
gtensor=[2.5 2.5 0; 2.5 2.5 0; 0 0 0]
sys = System(Nd_subcryst, dim, [SpinInfo(1, S=1/2,g=gtensor)], :SUN, units=Units.meV)
############## model set up ##############
J₁ = -0.33  #true value=-0.335*2
set_exchange!(sys, [J₁  -J₁   0.0;
                   -J₁   J₁   0.0;
                    0.0  0.0  0.0], Bond(3, 6, [0, 0, 0]))
# set_exchange!(sys, -0.5, Bond(3, 6, [0, 0, 0]))

J₂ = 0.285
set_exchange!(sys,[J₂   J₂ 0.0;
                   J₂   J₂ 0.0;
                   0.0  0.0 0.0], Bond(1, 7, [0, 0, 0]))

set_exchange!(sys,[J₂   J₂ 0.0;
                   J₂   J₂ 0.0;
                   0.0  0.0 0.0], Bond(7, 10, [0, 0, 0]))


field=[]
magnetization=[]
for B in 0:0.1:4
    set_external_field!(sys, (B/sqrt(2), -B/sqrt(2), 0))
    kT = 0.0
    randomize_spins!(sys)
    nsweeps = 10_000
    sampler = LocalSampler(kT=kT, propose=@mix_proposals 0.9 propose_uniform 0.1 propose_flip)
    for i in 1:nsweeps
       step!(sys, sampler)
    end
    sys_swt = reshape_geometry(sys, diagm(ones(Int, 3)))
    kT = 0.00
    randomize_spins!(sys)
    nsweeps = 100_000
    sampler = LocalSampler(kT=kT, propose=@mix_proposals 0.9 propose_uniform 0.1 propose_flip)
    for i in 1:nsweeps
        step!(sys_swt, sampler)
    end
    mx=0
    my=0
    mz=0
    for site in 1:length(sys_swt.dipoles)
        mx=mx+Sunny.magnetic_moment(sys_swt,site)[1]
        my=my+Sunny.magnetic_moment(sys_swt,site)[2]
        my=my+Sunny.magnetic_moment(sys_swt,site)[3]
    end
    mx=mx/length(sys_swt.dipoles)
    my=my/length(sys_swt.dipoles)
    mz=mz/length(sys_swt.dipoles)
    m=sqrt(mx^2+my^2+mz^2)
    append!(field,B)
    append!(magnetization,m)
end
Plots.plot(field,magnetization,line = (:steppre, :dot, :arrow, 0.5, 4, :red))