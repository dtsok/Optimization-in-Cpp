# Optimization-in-Cpp
Mathematical optimization project implementing derivative-free methods for minimization of the objective function.

### Task:
Minimize Lennard-Jones potential function $E_{LJ}$
$$E_{LJ} = 4\epsilon\sum_{i=1}^{N-1}\sum_{j=i+1}^{N}\Big[\big(\frac{\sigma}{r_{ij}}\big)^{12}-\big(\frac{\sigma}{r_{ij}}\big)^{6}\Big],$$
where $N$ is the number of the clusters, $\epsilon\ (=1)$ the *pair well depth*, $\sigma\ (=1)$ the *equilibrium pair separation* and $r_{ij}$ the Euclidean distance between $i$ and $j$ atom.

---

### Methods used:
1. Nelder-Mead Algorithm
2. Genetic Algorithm with binary representation
3. Genetic Algorithm with real representation
4. Particle Swarm Optimization

### Compile & Run
`make clean && make`

`./run <int>`

(not complete)
