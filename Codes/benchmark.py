import subprocess
import time
import matplotlib.pyplot as plt

def run_benchmark(files, num_threads, repetitions):
    results = []

    for filename in files:
        for num_thread in num_threads:
            for _ in range(repetitions):
                start_time = time.time()
                process = subprocess.Popen(
                    ["./bassin", filename],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env={"OMP_NUM_THREADS": str(num_thread)}
                )
                stdout, stderr = process.communicate()
                elapsed_time = (time.time() - start_time) * 1000 

                results.append({
                    "filename": filename,
                    "num_threads": num_thread,
                    "time": elapsed_time
                })
    
    return results

def plot_results(results):
    plt.figure(figsize=(10, 6))

    files = sorted(set(r["filename"] for r in results))
    num_threads = sorted(set(r["num_threads"] for r in results))

    for file in files:
        avg_times = []
        for num_thread in num_threads:
            times = [r["time"] for r in results if r["filename"] == file and r["num_threads"] == num_thread]
            avg_times.append(sum(times) / len(times))

        plt.plot(num_threads, avg_times, marker='o', label=file)

    plt.title("Benchmark OpenMP")
    plt.xlabel("Nombre de threads")
    plt.ylabel("Temps d'exécution moyen (ms)")
    plt.legend()
    plt.grid(True)
    plt.savefig("benchmark_openmp_graph.png")
    plt.show()

files = ["alpes.txt", "jeu_essai.txt", "grd_618360_6754408.txt"]
num_threads = [2, 4, 6, 8]
repetitions = 30

results = run_benchmark(files, num_threads, repetitions)
plot_results(results)
