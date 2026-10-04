// Version C: optimizada, cada pareja (i, j) con j > i se calcula una sola vez.
// Compilar: g++ -O2 -std=c++17 -fopenmp src/nbody_optimized.cpp -o nbody_optimized
// Ejecutar: ./nbody_optimized N steps threads schedule chunk

// INICIO BLOQUE COMUN
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include <omp.h>

const char* STUDENT = "Nombre1 Apellido1, Nombre2 Apellido2, Nombre3 Apellido3";

struct Body {
    double x, y;
    double vx, vy;
    double mass;
};

const double G = 1.0;          // constante gravitacional (unidades de simulacion)
const double DT = 0.01;        // paso de tiempo
const double EPS2 = 0.01;      // softening al cuadrado (evita dividir entre cero)
const unsigned SEED = 42;

struct Config {
    int n = 5000;
    int steps = 10;
    int threads = 1;
    std::string schedule = "static";
    int chunk = 0;             // 0 = chunk por defecto
};

// Argumentos en orden fijo: N steps threads schedule chunk
Config parse_args(int argc, char** argv) {
    Config c;
    if (argc > 1) c.n = std::atoi(argv[1]);
    if (argc > 2) c.steps = std::atoi(argv[2]);
    if (argc > 3) c.threads = std::atoi(argv[3]);
    if (argc > 4) c.schedule = argv[4];
    if (argc > 5) c.chunk = std::atoi(argv[5]);
    return c;
}

// Mismas condiciones iniciales en todas las versiones (seed = 42)
std::vector<Body> init_bodies(int n) {
    std::mt19937 gen(SEED);
    std::uniform_real_distribution<double> pos(0.0, 100.0);
    std::uniform_real_distribution<double> vel(-1.0, 1.0);
    std::uniform_real_distribution<double> mass(1.0, 10.0);
    std::vector<Body> b(n);
    for (int i = 0; i < n; i++) {
        b[i].x = pos(gen);
        b[i].y = pos(gen);
        b[i].vx = vel(gen);
        b[i].vy = vel(gen);
        b[i].mass = mass(gen);
    }
    return b;
}

// Aceleracion = fuerza / masa; luego se actualiza velocidad y posicion
void update_bodies(std::vector<Body>& b, const std::vector<double>& fx, const std::vector<double>& fy) {
    for (size_t i = 0; i < b.size(); i++) {
        b[i].vx += fx[i] / b[i].mass * DT;
        b[i].vy += fy[i] / b[i].mass * DT;
        b[i].x += b[i].vx * DT;
        b[i].y += b[i].vy * DT;
    }
}

// Checksum: suma de las posiciones finales
double checksum(const std::vector<Body>& b) {
    double s = 0.0;
    for (const Body& p : b) s += p.x + p.y;
    return s;
}

void print_summary(const char* exe, const char* mode, const Config& c, int threads_used,
                   double elapsed, double sum) {
    std::printf("Student: %s\n", STUDENT);
    std::printf("Executable: %s\n", exe);
    std::printf("Mode: %s\n", mode);
    std::printf("N: %d | Steps: %d | Seed: %u\n", c.n, c.steps, SEED);
    std::printf("Threads requested: %d | Threads used: %d\n", c.threads, threads_used);
    if (c.chunk > 0)
        std::printf("Schedule: %s | Chunk: %d\n", c.schedule.c_str(), c.chunk);
    else
        std::printf("Schedule: %s | Chunk: default\n", c.schedule.c_str());
    std::printf("Elapsed: %.3f s\n", elapsed);
    std::printf("Checksum: %.6f\n", sum);
    std::printf("Processors reported by OpenMP: %d\n", omp_get_num_procs());
    std::printf("Max threads reported by OpenMP: %d\n", omp_get_max_threads());
}
// FIN BLOQUE COMUN

int main(int argc, char** argv) {
    Config cfg = parse_args(argc, argv);

    // Threads, schedule y chunk vienen de los argumentos; el pragma usa schedule(runtime)
    omp_sched_t kind = omp_sched_static;
    if (cfg.schedule == "dynamic") kind = omp_sched_dynamic;
    else if (cfg.schedule == "guided") kind = omp_sched_guided;
    omp_set_num_threads(cfg.threads);
    omp_set_schedule(kind, cfg.chunk);   // chunk = 0 usa el valor por defecto

    // Threads realmente usados por una region paralela
    int threads_used = 1;
    #pragma omp parallel
    {
        #pragma omp single
        threads_used = omp_get_num_threads();
    }

    std::vector<Body> bodies = init_bodies(cfg.n);
    std::vector<double> fx(cfg.n), fy(cfg.n);
    // Acumuladores privados: un arreglo de fuerzas por thread (thread t usa [t * n, t * n + n))
    std::vector<double> fxt((size_t)threads_used * cfg.n), fyt((size_t)threads_used * cfg.n);

    double t0 = omp_get_wtime();
    for (int s = 0; s < cfg.steps; s++) {
        #pragma omp parallel
        {
            double* px = &fxt[(size_t)omp_get_thread_num() * cfg.n];
            double* py = &fyt[(size_t)omp_get_thread_num() * cfg.n];
            for (int i = 0; i < cfg.n; i++) px[i] = py[i] = 0.0;

            // Fase 1a: cada pareja una vez; se suma a i y se resta a j (tercera ley de Newton).
            // Varios threads tocan el mismo j, por eso cada uno escribe en su propio arreglo.
            #pragma omp for schedule(runtime)
            for (int i = 0; i < cfg.n; i++) {
                double sx = 0.0, sy = 0.0;
                for (int j = i + 1; j < cfg.n; j++) {
                    double dx = bodies[j].x - bodies[i].x;
                    double dy = bodies[j].y - bodies[i].y;
                    double r2 = dx * dx + dy * dy + EPS2;
                    double f = G * bodies[i].mass * bodies[j].mass / (r2 * std::sqrt(r2));
                    sx += f * dx;
                    sy += f * dy;
                    px[j] -= f * dx;
                    py[j] -= f * dy;
                }
                px[i] += sx;
                py[i] += sy;
            }
            // Barrera implicita: todos los threads terminaron sus parejas

            // Fase 1b: sumar los arreglos de todos los threads; cada i lo escribe un solo thread
            #pragma omp for schedule(static)
            for (int i = 0; i < cfg.n; i++) {
                double sx = 0.0, sy = 0.0;
                for (int t = 0; t < threads_used; t++) {
                    sx += fxt[(size_t)t * cfg.n + i];
                    sy += fyt[(size_t)t * cfg.n + i];
                }
                fx[i] = sx;
                fy[i] = sy;
            }
        }
        // Fase 2: actualizar velocidades y posiciones
        update_bodies(bodies, fx, fy);
    }
    double elapsed = omp_get_wtime() - t0;

    print_summary(argv[0], "optimized", cfg, threads_used, elapsed, checksum(bodies));
    return 0;
}
