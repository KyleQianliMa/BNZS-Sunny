using Sunny
using LinearAlgebra
using GLMakie, Plots

############## load cif file ##############
cif_file = joinpath(@__DIR__, "BNZS_SSL.cif")
cryst = Crystal(cif_file, symprec=1e-3)
Nd_subcryst = subcrystal(cryst, "Nd1");
# view_crystal(Nd_subcryst, 7.0)

############## symmetry analysis ##############
# print_symmetry_table(Nd_subcryst, 5.0)
dim = (1, 1, 1)
gtensor=[2.6 2.6 0; 2.6 2.6 0; 0 0 0]
sys = System(Nd_subcryst, dim, [SpinInfo(1, S=1/2,g=gtensor)], :SUN, units=Units.meV)

# print_bond(Nd_subcryst, Bond(3, 17, [0, 0, 0]); b_ref=Bond(1,19,[0,0,0]))

############## model set up ##############
J₁ = -0.33  #true value=-0.335*2
set_exchange!(sys, [J₁  -J₁   0.0;
                   -J₁   J₁   0.0;
                    0.0  0.0  0.0], Bond(3, 6, [0, 0, 0]))

J₂ = 0.285
set_exchange!(sys,[J₂   J₂ 0.0;
                   J₂   J₂ 0.0;
                   0.0  0.0 0.0], Bond(1, 7, [0, 0, 0]))

set_exchange!(sys,[J₂   J₂ 0.0;
                   J₂   J₂ 0.0;
                   0.0  0.0 0.0], Bond(7, 10, [0, 0, 0]))

B = 2
set_external_field!(sys, (B/sqrt(2), -B/sqrt(2), 0))

############## use Monte Carlo to find the ground state magnetic ordering ##############
kT = 0.0
randomize_spins!(sys)
nsweeps = 10_000
sampler = LocalSampler(kT=kT, propose=@mix_proposals 0.5 propose_uniform 0.5 propose_flip)
e=[]
n=[]
for i in 1:nsweeps
    step!(sys, sampler)
    append!(e,energy(sys))
    append!(n,i)
end
Plots.plot(n,e,line = (:steppre, :dot, :arrow, 0.5, 4, :blue))

################ analyze the ground state ##############
print_wrapped_intensities(sys)
suggest_magnetic_supercell([[0, 0, 0]], sys.latsize) #copy printed intensity wavevectors
sys_swt = reshape_geometry(sys, diagm(ones(Int, 3)))

# suggest_magnetic_supercell([[0, 0, 0],[1/2, 1/2, 0],[1/2, 0, 0],[0, 1/2, 0]], sys.latsize)
# sys_swt = reshape_geometry(sys,  [2 0 0; 0 2 0; 0 0 1])

################ further optimize the ground state for a spin-wave system ##############
kT = 0.0
randomize_spins!(sys)
nsweeps = 100_000
sampler = LocalSampler(kT=kT, propose=@mix_proposals 0.3 propose_uniform 0.7 propose_flip)
e=[]
n=[]
for i in 1:nsweeps
    step!(sys_swt, sampler)
    append!(e,energy(sys_swt))
    append!(n,i)
end
Plots.plot(n,e,line = (:steppre, :dot, :arrow, 0.5, 4, :red))

plot_spins(sys_swt,arrowlength=1, linewidth=0.5, arrowsize=1.0)

################ spin-wave calculations ##############
swt = SpinWaveTheory(sys_swt)
qvals = 0.0:0.005:1.0
qs = [[q, q, 0] for q in qvals]
lenq = length(qs)
numband = length(sys_swt.dipoles)
disp = dispersion(swt, qs)
energies = 0.0:0.01:2.0
INS_intensities = Sunny.intensities(swt, qs, energies, 0.02)
# Plots.plot(qvals, disp, linecolor="blue", label="", xlims=(0, 1), ylims=(0, 2), xlabel="(q, q, 0)", ylabel="E (meV)")
Plots.heatmap(qvals, energies, INS_intensities', xlims=(0, 1), ylims=(0, 2), xlabel="(q, q, 0)", ylabel="E (meV)", title=B)