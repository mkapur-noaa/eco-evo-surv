# Load required package
# using Pkg
#Pkg.add("Distributions")
#Pkg.add("StatsBase")
#Pkg.add("PDMats")
#Pkg.add("NeutralLandscapes")
#Pkg.add(url = "https://github.com/jkbest2/TweedieDistributions.jl")
#Pkg.add(url = "https://github.com/jkbest2/FisherySim.jl")
using FisherySim
using Distributions
using StatsBase
using Test
using Random
using PDMats
using LinearAlgebra
using NeutralLandscapes

## Make reproducible
Random.seed!(1234)

## -----------------------------------------------------------------------------
## Construct a GriddedFisheryDomain
origin = (0.0, 0.0)
antipode = (100.0, 100.0)
n = (50, 50)
Ω = GriddedFisheryDomain(origin, antipode, n)
T = 10 # Number of years

include("julia/jkbest2/test-fisherydomain.jl")

## -----------------------------------------------------------------------------
## Construct covariance kernels and covariance matrices over fishery domains
σ² = 5.0
ρ = 40.0

expkern = ExpCov(σ², ρ)
matkern = Matérn32Cov(σ², ρ)
ar1kern = AR1(0.1, 0.5)

expΣ = cov(expkern, Ω)
matΣ = cov(matkern, Ω)
ar1Σ = cov(ar1kern, T)

include("test-covkernels.jl")

## -----------------------------------------------------------------------------
## Define distributions over the domain
lognorm = DomainDistribution(LogNormal(-0.2^2 / 2, 0.2), Ω)
mvlognorm = DomainDistribution(MvLogNormal, ones(length(Ω)), matΣ, Ω)
mvnorm = DomainDistribution(MvNormal(zeros(length(Ω)), matΣ), Ω)
# matlognorm = DomainDistribution(MatrixLogNormal, 0.2ones(length(Ω), T), matΣ, ar1Σ, Ω)

ln_real = rand(lognorm)
mvln_real = rand(mvlognorm)
mvn_real = rand(mvnorm)
# matlognorm_real = rand(matlognorm)

dd1 = DomainDistribution(MidpointDisplacement(), Ω)
cdd1 = ClassifiedDomainDistribution(MidpointDisplacement(), Ω, 3)
cdd2 = ClassifiedDomainDistribution(MidpointDisplacement(), Ω, [0.1, 0.4])
ddist_vec = [mvnorm, dd1]
mdd1 = MultiDomainDistribution(ddist_vec)
bdd1 = BlendedDomainDistribution(mdd1, [0.3, 0.5])
bdd2 = BlendedDomainDistribution(ddist_vec, Ω, [0.5, 0.9])
bdd3 = BlendedDomainDistribution([mvlognorm, cdd1], Ω, [0.3, 2.0])

include("julia/jkbest2/test-domaindistribution.jl")

## -----------------------------------------------------------------------------
## Habitats
hab1 = Habitat(rand(dd1))
hab2 = Habitat(rand.(ddist_vec))
include("julia/jkbest2/test-habitat.jl")

## -----------------------------------------------------------------------------
## Generate bathymetry
μv = zeros(length(Ω))
# μm = zeros(size(Ω)...)
for idx in eachindex(Ω)
    μv[idx] = -1 / 2500 * first(Ω[idx]) * (first(Ω[idx]) - 100)
end

# hab_model = DomainDistribution(MvLogNormal, μv, matΣ, Ω)
hab_model = DomainDistribution(MvNormal(μv, matΣ), Ω)
hab = rand(hab_model)

hab2 = Habitat(rand(mdd1))

## -----------------------------------------------------------------------------
## Construct a movement model
h1pref(h) = exp(-h ^ 2 / 10)
h2pref(h) = 1
hab2pref = HabitatPreference(h1pref, h2pref)
dist = MovementRate(d -> exp(-d ^ 2 / 100))
move = MovementModel(hab2, hab2pref, dist, Ω)

# Find the equilibrium distribution using `eigvecs` to compare to to approximate
# method.
ews = eigvecs(move.M)
eqdist_eigvecs = real.(reshape(ews[:, end], size(Ω)...))
eqdist_eigvecs ./= sum(eqdist_eigvecs)

eqdist_ap0 = eqdist(move)
eqdist_ap = eqdist(move, 100.0)

include("julia/jkbest2/test-movement.jl")

## -----------------------------------------------------------------------------
## Define Schaefer and PellaTomlinson dynamics models
r = 0.05
K = 100.0

schaef = Schaefer(r, K)
P1 = schaef(eqdist_ap)
Phalf = PopState(P1.P ./ 2)
Phalf1 = schaef(Phalf)

pt = PellaTomlinson(r, K, 3.39)
P1_pt = pt(eqdist_ap)
P1half_pt = pt(Phalf1)

unispt = StochasticProduction(PellaTomlinson(r, K, 3.39),
                              LogNormal(-0.2^2 / 2, 0.2))
multispt = StochasticProduction(
    PellaTomlinson(r, K, 3.39),
    mvlognorm)
P1_uspt = unispt(P1)
P1_mspt = multispt(P1)

include("julia/jkbest2/test-pop_dynamics.jl")

## -----------------------------------------------------------------------------
## Define vessels and fleets
target_rand = RandomTargeting()
survey_stations = vec(LinearIndices(n)[2:2:48, 2:8:48])
target_fixed = FixedTargeting(survey_stations)
target_pref = PreferentialTargeting(10 .* eqdist_ap.P, Ω)

fixed_t = target(target_fixed, Ω)
reset!(target_fixed)
rand_t = target(target_rand, Ω)
reset!(target_rand)
pref_t = target(target_pref, Ω)
reset!(target_pref)

q_const = Catchability(0.2)
q_diff = Catchability(lognorm, q_const)
q_spat = Catchability(mvlognorm, q_const)
# q_sptemp = Catchability(matlognorm, q_const)
q_hab = HabitatCatchability(hab2, 0.1, h1 -> 1.0, h2 -> 1.2 * h2)

## These have to be fairly high to get many non-zero catches.
## These give ~10% zeros. These (especially ϕ?) may be good for testing,
## as they result in completely fishing out some cells, testing that check.
ξ = 1.9
ϕ = 1.9

v1 = Vessel(target_fixed, q_const, ξ, ϕ)
v2 = Vessel(target_rand, q_diff, ξ, ϕ)
v3 = Vessel(target_pref, q_spat, ξ, ϕ)
v4 = Vessel(target_pref, q_hab, ξ, ϕ)

P = PopState(ones(size(Ω)...) / length(Ω))
c1 = fish!(P, v1, Ω)
reset!(v1.target)

fleet = Fleet([v1, v2, v3, v4], [length(v1.target), 100, 1000, 1000])

c2 = fish!(P, fleet, Ω)
total_catch = sum(getfield.(c2, :catch_biomass))

include("julia/jkbest2/test-vessels.jl")

## -----------------------------------------------------------------------------
## Fish a population using above definitions
P1 = eqdist_ap
Pvec, Cvec = simulate(P1, fleet, move, schaef, Ω, T)

Psums = sum.(Pvec)

function getcatch(CV::Vector{C}) where {C<:Catch}
    timevec = getfield.(CV, :time)
    catchvec = getfield.(CV, :catch_biomass)
    catchtot = zeros(length(unique(timevec)))
    for (t, c) in zip(timevec, catchvec)
        catchtot[t] += c
    end
    catchtot
end

Csums = getcatch(Cvec)

p2 = begin
         rempop = Psums[1] - Csums[1]
         rempop + rempop * r * (1 - rempop / K)
     end

include("julia/jkbest2/test-simulate.jl")

