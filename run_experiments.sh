#!/bin/bash
# Compila las tres versiones y corre cada configuracion 3 veces.
# Los resultados se guardan en results.csv.
N=5000
STEPS=10
RUNS=3

g++ -O2 -std=c++17 -fopenmp src/nbody_seq.cpp -o nbody_seq || exit 1
g++ -O2 -std=c++17 -fopenmp src/nbody_parallel.cpp -o nbody_parallel || exit 1
g++ -O2 -std=c++17 -fopenmp src/nbody_optimized.cpp -o nbody_optimized || exit 1

echo "version,threads,schedule,chunk,run,elapsed,checksum" > results.csv

# Ejecuta un binario y agrega una fila al CSV
run() {
    out=$(./$1 $N $STEPS $3 $4 $5)
    elapsed=$(echo "$out" | grep "Elapsed:" | awk '{print $2}')
    sum=$(echo "$out" | grep "Checksum:" | awk '{print $2}')
    echo "$2,$3,$4,$5,$r,$elapsed,$sum" >> results.csv
    echo "$2 threads=$3 $4 chunk=$5 run=$r -> $elapsed s"
}

# Configuraciones de la seccion 5: threads schedule chunk (0 = default)
CONFIGS=(
    "2 static 0"
    "4 static 0"
    "8 static 0"
    "8 static 8"
    "8 static 64"
    "8 dynamic 8"
    "8 dynamic 64"
    "8 guided 0"
)

for r in $(seq 1 $RUNS); do
    run nbody_seq sequential 1 - 0
    for c in "${CONFIGS[@]}"; do
        run nbody_parallel parallel $c
        run nbody_optimized optimized $c
    done
done
