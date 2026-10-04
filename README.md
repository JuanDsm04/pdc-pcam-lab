# Laboratorio PCAM + OpenMP: Simulación N-Body 2D

## Compilar

```bash
g++ -O2 -std=c++17 -fopenmp src/nbody_seq.cpp -o nbody_seq
g++ -O2 -std=c++17 -fopenmp src/nbody_parallel.cpp -o nbody_parallel
g++ -O2 -std=c++17 -fopenmp src/nbody_optimized.cpp -o nbody_optimized
```

## Ejecutar

Los tres programas reciben los mismos argumentos, en este orden:

```bash
./nbody_parallel N steps threads schedule chunk
```

- `schedule`: `static`, `dynamic` o `guided`.
- `chunk`: tamaño de chunk; `0` usa el valor por defecto.
- La versión secuencial ignora `threads`, `schedule` y `chunk`.

Ejemplo: `./nbody_parallel 5000 10 8 dynamic 64`

## Experimentos

```bash
bash run_experiments.sh
```

Compila las tres versiones, corre cada configuración de la sección 5 tres veces y guarda los resultados en `results.csv`.

## Checksum de referencia

Versión secuencial, N = 5000, 10 pasos, seed = 42 (g++ en Linux):

```
Checksum: 497859.295223
```

Las distribuciones aleatorias de C++ pueden variar entre compiladores, así que las versiones siempre se comparan en la misma máquina.
