include("scr/JuPCHS.jl")
using .JuPCHS

# Parameters
l = 15.91e-3    # H
c = 50e-6       # F
r = 25          # ohm
v0 = 24         # v

# Port-Hamiltonian formulation
Q = [1/l 0; 0 1/c]
Js = [0 -1; 1 0]
Rs = [0 0; 0 1/r]
Gs = [v0; 0]
Hs(x) = 0.5*(x'*(Q*x))
dHs(x) = Q*x
buck = BuildPCHS(nx=2,nu=1, dt=1/45e3, 
                 R= Rs, J=Js, H=Hs, dH=dHs,G=Gs,
                 name="buck",
                 xlabel=["li","cv"],
                 ulabel=["m"])

# Simulation
u_input(x,t) = [0.5*Heaviside(t-0.02)]
results = Simulate(PCHS=buck,nt=4000,xini=[0;0],u=u_input,method="Euler")
PlotResults(PCHS=buck,data=results,graph=["x","H"],xscale=[1/l;1/c])
            