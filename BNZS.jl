using Sunny
using LinearAlgebra
using GLMakie, Plots

############## load cif file ##############
cif_file = joinpath(@__DIR__, "BNZS_SSL.cif")
cryst = Crystal(cif_file, symprec=1e-3)
Nd_subcryst = subcrystal(cryst, "Nd1")
# view_crystal(Nd_subcryst, 7.0)

############## symmetry analysis ##############
print_symmetry_table(Nd_subcryst, 5.0)
dim = (1, 1, 1)
sys = System(Nd_subcryst, dim, [SpinInfo(1, S=1/2)], :SUN)

############## model set up ##############
J₁ = -0.66
set_exchange!(sys, [J₁/2  -J₁/2   0.0;
                   -J₁/2   J₁/2   0.0;
                    0.0  0.0  0.0], Bond(3, 6, [0, 0, 0]))

J₂ = 0.59
set_exchange!(sys,[J₂/2   J₂/2 0.0;
                   J₂/2   J₂/2 0.0;
                   0.0  0.0 0.0], Bond(1, 7, [0, 0, 0]))

set_exchange!(sys,[J₂/2   J₂/2 0.0;
                   J₂/2   J₂/2 0.0;
                   0.0  0.0 0.0], Bond(7, 10, [0, 0, 0]))
B = 1
set_external_field!(sys, (sqrt(B/2), -sqrt(B/2), 0))

############## use Monte Carlo to find the ground state magnetic ordering ##############
kT = 0
randomize_spins!(sys)
nsweeps = 100_000
sampler = LocalSampler(kT=kT, propose=propose_uniform)
for i in 1:nsweeps
    step!(sys, sampler)
end

plot_spins(sys,arrowlength=1, linewidth=0.5, arrowsize=1.0)

################ analyze the ground state ##############
print_wrapped_intensities(sys)
suggest_magnetic_supercell([[0, 0, 0]], sys.latsize)
sys_swt = reshape_geometry(sys, diagm(ones(Int, 3)))

################ further optimize the ground state for a spin-wave system ##############
kT = 0
randomize_spins!(sys)
nsweeps = 50_000
sampler = LocalSampler(kT=kT, propose=propose_uniform)
for i in 1:nsweeps
    step!(sys_swt, sampler)
end

################ spin-wave calculations ##############
swt = SpinWaveTheory(sys_swt)
qvals = 0.0:0.005:1.0
qs = [[q, q, 0] for q in qvals]
lenq = length(qs)
numband = length(sys_swt.dipoles)
disp = dispersion(swt, qs)
energies = 0.0:0.01:2.0
INS_intensities = Sunny.intensities(swt, qs, energies, 0.05)
Plots.plot(qvals, disp, linecolor="blue", label="", xlims=(0, 1), ylims=(0, 2), xlabel="(q, q, 0)", ylabel="E (meV)")
Plots.heatmap(qvals, energies, INS_intensities', xlims=(0, 1), ylims=(0, 2), xlabel="(q, q, 0)", ylabel="E (meV)")