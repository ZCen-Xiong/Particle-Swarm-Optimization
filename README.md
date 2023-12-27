# Particle Swarm Algorithm (Without External Dependencies) 
Built-in matrix operation function, does not depend on other function libraries, the optimization process is as follows:

$$
\begin{algorithm}
\caption{Particle Swarm Optimization (PSO)}
\begin{algorithmic}[1]
\State Initialize population of particles with random positions and velocities
\State Initialize global best position vector: $P_{\text{global}}$
\While{stopping criteria are not met}
    \For{each particle}
        \State Update particle’s best position vector: $p_i$
        \State Update global best position vector: $P_{\text{global}}$
        \For{each dimension}
            \State Update particle’s velocity: $v_i(t + 1) = w * v_i(t) + c1 * r1 * (p_i(t) - x_i(t)) + c2 * r2 * (P_{\text{global}}(t) - x_i(t))$
            \State Update particle’s position: $x_i(t + 1) = x_i(t) + v_i(t + 1)$
        \EndFor
    \EndFor
\EndWhile
\State \Return global best position: $P_{\text{global}}$
\end{algorithmic}
\end{algorithm}
$$
