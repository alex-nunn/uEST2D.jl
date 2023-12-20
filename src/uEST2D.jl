module uEST2D

using DifferentialEquations

export TF_ϕ0, TF_L, TF_lengths, TF_potentials, TF_GAP
export ϕ_electrode, E_electrode!, trap_potential, trap_E!
export TF_trap_potential, TF_trap_E!
export trap_signal
export ion_timescale, ion_trajectory, poincare_section

# ------------------------------------------------------------------------------
#  Fixed constants
# ------------------------------------------------------------------------------
# Characteristic scales
const TF_ϕ0 = 100  # Volts
const TF_L = 60    # μm

# Thermo Fisher trap configuration
# -- electrode lengths
const TF_lengths_dim = [
    188,  # μm - Electrode 1 - center Einzel lens
    454,  # μm - Electrode 2 - pickup electrode
    140,  # μm - Electrode 3 - reflector Einzel lens
    45,   # μm - Electrode 4 - reflector
    45,   # μm - Electrode 5 - ...
    21,   # μm - Electrode 6 - ...
    124   # μm - Electrode 7 - outer reflector
]

# -- electrode potentials
const TF_potentials_dim = [
    34.03,    # Volts - Electrode 1 - center Einzel lens
    0.00,     # Volts - Electrode 2 - pickup electrode
    -143.01,  # Volts - Electrode 3 - reflector Einzel lens
    10.09,    # Volts - Electrode 4 - reflector
    127.06,   # Volts - Electrode 5 - ...
    146.80,   # Volts - Electrode 6 - ...
    178.15    # Volts - Electrode 7 - outer reflector
]

const TF_GAP = 0.0     # μm  -- gap between adjacent electrodes

# dimensionless configuration values
const TF_lengths = TF_lengths_dim ./ TF_L
const TF_potentials = TF_potentials_dim ./ TF_ϕ0


# ------------------------------------------------------------------------------
# Electrostatic field
# ------------------------------------------------------------------------------
"""
    ϕ_electrode(y, z, a)

Potential due to a 2D electrode of length `a` in yz-plane positioned at 
y = 0, z ∈ (0, a) where the 2D trap has planes at y = 0, 1.
"""
ϕ_electrode(y, z, a) = 0.0 <= y <= 1.0 ? (
    atan(cot(π * y / 2) * tanh(π * (a - z) / 2)) 
    + atan(cot(π * y / 2) * tanh(π * z / 2)) ) / π : 0.0


"""
    E_electrode!(out, y, z, a, ϕ0, s)

Add 2D electric field due to an electrode of length `a` in yz-plane at
y = 0, z ∈ (0, a) with potential `ϕ0` to `out`. The `s` factor controls the sign 
of the electric field in the `y` direction, which is useful for modelling
electrodes on the top face of the μEST.
"""
E_electrode!(out, y, z, a, ϕ0, s) = begin
    coshπaz = cosh(π * (a - z))
    coshπz = cosh(π * z)
    cosπy = cos(π * y) 
    C1 = cosπy - coshπaz
    C2 = cosπy - coshπz

    out[1] += -ϕ0 * s * 0.5(sinh(π * (a - z)) / C1 + sinh(π * z) / C2)
    out[2] += -ϕ0 * 0.5(coshπaz -  coshπz) * sin(π * y) / (C1 * C2)
end


"""
    trap_potential(y, z, electrode_lengths, electrode_pots, gap, axial_d)

Return electrostatic potential within 2D μEST at (y, z)

# Arguments
- `y`: vertical position in trap in dimensionless units, y ∈ (-1/2, 1/2)
- `z`: horizontal postion in trap in dimensionless units z ∈ (-∞, ∞)
- `electrode_lengths::Vector`: lengths of electrodes on top and bottom faces
- `electrode_pots::Vector`: potentials of electrodes on top and bottom faces
- `gap`: gap between adjacent electrodes
- `axial_d`: displacement of top face relative to bottom in z-direction
"""
function trap_potential(y, z, electrode_lengths, electrode_pots, gap, axial_d=0.0)
    ϕ1 = electrode_pots[begin]
    l1 = electrode_lengths[begin]
    
    value = ϕ1 * (
        ϕ_electrode(y + 0.5, z + 0.5l1, l1) 
        + ϕ_electrode(0.5 - y, z + 0.5l1 - axial_d, l1)
    )

    z_pos = 0.5l1 + gap
    for (ln, ϕn) ∈ zip(electrode_lengths[2:end], electrode_pots[2:end])
        value += ϕn * (
            ϕ_electrode(y + 0.5, z - z_pos, ln)
            + ϕ_electrode(0.5 - y, z - z_pos - axial_d, ln)
            + ϕ_electrode(y + 0.5, z + z_pos + ln, ln)
            + ϕ_electrode(0.5 - y, z + z_pos + ln - axial_d, ln)
        )
        z_pos += ln + gap
    end

    return value
end


"""
    trap_E!(out, pt, electrode_lengths, electrode_pots, gap, axial_d)

Add electrostatic field within 2D μEST at position `pt` to `out` vector

# Arguments
- `out`: electric field output
- `pt`: position within trap `pt = [y, z]` where y ∈ (-1/2, 1/2), z ∈ (-∞, ∞)
- `electrode_lengths::Vector`: lengths of electrodes on top and bottom faces
- `electrode_pots::Vector`: potentials of electrodes on top and bottom faces
- `gap`: gap between adjacent electrodes
- `axial_d`: displacement of top face relative to bottom in z-direction
"""
function trap_E!(out, pt, electrode_lengths, electrode_pots, gap, axial_d=0.0)
    ϕ1 = electrode_pots[begin]
    l1 = electrode_lengths[begin]
    
    y, z = pt
    out[:] .= 0.0

    E_electrode!(out, y + 0.5, z + 0.5l1, l1, ϕ1, 1)
    E_electrode!(out, 0.5 - y, z + 0.5l1 - axial_d, l1, ϕ1, -1)

    z_pos = 0.5l1 + gap

    for i ∈ eachindex(electrode_lengths)[2:end]
        ln = electrode_lengths[i]
        ϕn = electrode_pots[i]
        
        E_electrode!(out, y + 0.5, z - z_pos, ln, ϕn, 1)
        E_electrode!(out, 0.5 - y, z - z_pos - axial_d, ln, ϕn, -1)
        E_electrode!(out, y + 0.5, z + z_pos + ln, ln, ϕn, 1)
        E_electrode!(out, 0.5 - y, z + z_pos + ln - axial_d, ln, ϕn, -1)

        z_pos += ln + gap
    end
end


"""
    TF_trap_potential(y, z)

Return electrostatic potential at (y, z) within Thermo Fisher μEST

See also [`trap_potential`](@ref)
"""
function TF_trap_potential(y, z, axial_d=0.0)
    return trap_potential(y, z, TF_lengths, TF_potentials, TF_GAP, axial_d)
end


"""
    TF_trap_E!(out, pt)

Add electrostatic field at `pt` to `out` for Thermo Fisher μEST

See also [`trap_E!`](@ref)
"""
function TF_trap_E!(out, pt, axial_d=0.0)
    trap_E!(out, pt, TF_lengths, TF_potentials, TF_GAP, axial_d)
end


# ------------------------------------------------------------------------------
# Pickup electrode signal
# ------------------------------------------------------------------------------
"""
    trap_signal(y, z, electrode_lengths, gap, axial_d)

Return signals for an ion at position (y, z). The components of the returned 
vector are the accumulated charge on the pickup electrodes as,

    [Q_bottom_left, Q_bottom_right, Q_top_left, Q_top_right]
"""
function trap_signal(y, z, axial_d=0.0, electrode_lengths=TF_lengths, gap=TF_GAP)
    l1 = electrode_lengths[begin]
    z_pos = 0.5l1 + gap

    l2 = electrode_lengths[2]

    return [
        ϕ_electrode(y + 0.5, z + z_pos + l2, l2),           # bottom-left
        ϕ_electrode(y + 0.5, z - z_pos, l2),                # bottom-right
        ϕ_electrode(0.5 - y, z + z_pos + l2 - axial_d, l2), # top-left
        ϕ_electrode(0.5 - y, z - z_pos - axial_d, l2)       # top-right
    ]
end

# ------------------------------------------------------------------------------
# Ion trajectory
# ------------------------------------------------------------------------------
"""
    ion_timescale(mz_thompson, trap_height)

Return characteristic timescale for an ion trajectory in nanoseconds
    
# Arguments
- `mz_thompson`: mass to charge ratio of ion in units of Thompsons
- `trap_height`: height of the 2D μEST in units of micrometers
"""
function ion_timescale(mz_thompson, trap_height=TF_L)
    mt_dalton = 1.660
    qt_electron = 1.602

    T0 = trap_height * sqrt(mz_thompson * mt_dalton / qt_electron) * 1e-2
    return T0
end


"""
    ion_trajectory(E0, y0, α0, T, trap_ϕ, trap_E)

Return trajectory for an ion with energy `E0`, initial position `(y, 0)` and
insertion angle `α0`.

Trajectories are always started at the midplane of the 2D μEST.

# Arguments
- `E0`: dimensionless energy of the ion (determines initial velocity)
- `y0`: initial vertical position of ion in the midplane
- `α0`: angle of initial velocity in radians from horizontal 
        (α0 = 0 --> aligned with axis of trap)
- `T`: length of time to simulate trajectory over (in dimensionless units)
- `trap_ϕ`: electrostatic potential of trap with call signature `trap_ϕ(y, z)`
- `trap_E`: electrostatic field of trap with call signature `trap_E!(out, pt)`

See also [`poincare_section`](@ref)
"""
function ion_trajectory(
        E0, y0, α0, T, trap_ϕ=TF_trap_potential, trap_E=TF_trap_E!; solve_kws...
        )
    u0 = [y0, 0.0]
    vel0 = sqrt(2 * (E0 - trap_ϕ(u0...)))
    du0 = vel0 .* [sin(α0), cos(α0)]

    # Setup ODE system
    tspan = (0, T)
    f! = (ddu, du, u, p, t) -> trap_E(ddu, u)

    prob = SecondOrderODEProblem(f!, du0, u0, tspan)
    sol = solve(
        prob, DPRKN6();
        abstol=1e-8,
        reltol=1e-8,
        maxiters=1e8,
        solve_kws...
    )
    
    return sol
end


"""
    poincare_section(E0, y0, α0, T, z0, trap_ϕ, trap_E!)

Return a poincare section for 2D μEST at plane `z = z0`

# Arguments
- `E0`: dimensionless energy of the ion (determines initial velocity)
- `y0`: initial vertical position of ion in the midplane
- `α0`: angle of initial velocity in radians from horizontal 
        (α0 = 0 --> aligned with axis of trap)
- `T`: length of time to simulate trajectory over (in dimensionless units)
- `z0`: position of Poincare section
- `trap_ϕ`: electrostatic potential of trap with call signature `trap_ϕ(y, z)`
- `trap_E`: electrostatic field of trap with call signature `trap_E!(out, pt)`

See also [`ion_trajectory`](@ref)
"""
function poincare_section(
        E0, y0, α0, T, z0=0.0, trap_ϕ=TF_trap_potential, trap_E=TF_trap_E!;
        abstol=1e-8, reltol=1e-8
        )
    # Setup initial conditions
    u0 = [y0, z0]
    vel0 = sqrt(2 * (E0 - trap_ϕ(u0...)))
    du0 = vel0 .* [sin(α0), cos(α0)]

    # Setup ODE system
    tspan = (0, T)
    f! = (ddu, du, u, p, t) -> trap_E(ddu, u)

    αs = []
    ys = []
    ts = []

    # Setup callbacks
    poincare_cb = ContinuousCallback(
        (u, t, integrator) -> u[4] - z0,  # condition: z == 0
        (integrator) -> begin  # affect : store (y, α, t)
            u = integrator.u
            push!(ys, u[3])
            push!(αs, atan(u[1] / u[2]))
            push!(ts, integrator.t)
        end
    )
    bounds_cb = DiscreteCallback(
        # condition: out of bounds?
        (u, t, integrator) -> !(abs(u[3]) < 1 && abs(u[4]) < 17),
        # affect : stop trajectory
        terminate!
    )
    callbacks = CallbackSet(poincare_cb, bounds_cb)
    
    prob = SecondOrderODEProblem(f!, du0, u0, tspan, callback=callbacks)
    sol = solve(
        prob, DPRKN6(),
        save_on=false,
        save_start=true,
        abstol=abstol,
        reltol=reltol,
        maxiters=1e8
    )
    
    u_final = sol.u[end]
    E_final = 0.5 * (u_final[1] ^2 + u_final[2]^2) + trap_ϕ(u_final[3:4]...)
    E_rel_error = E_final / E0 - 1

    is_confined = sol.retcode == ReturnCode.Success
    return (
        y = ys, 
        α = αs, 
        t = ts, 
        E_error = E_rel_error, 
        is_confined = is_confined,
        E0 = E0,
        y0 = y0,
        α0 = α0,
        T = T
    )
end

end # end module
