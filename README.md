# Particle Swarm Optimization (Without External Dependencies) 
Built-in matrix operation function, does not depend on other function libraries, the optimization process is as follows:

## Particle Swarm Optimization (PSO)

1. Initialize population of particles with random positions and velocities
2. Initialize global best position vector: P_global
3. Repeat until stopping criteria are met:
    1. For each particle:
        1. Update particle’s best position vector: p_i
        2. Update global best position vector: P_global
        3. For each dimension:
            1. Update particle’s velocity:
                ```
                v_i(t + 1) = w * v_i(t) + c1 * r1 * (p_i(t) - x_i(t)) + c2 * r2 * (P_global(t) - x_i(t))
                ```
            2. Update particle’s position:
                ```
                x_i(t + 1) = x_i(t) + v_i(t + 1)
                ```
4. Return global best position: P_global
