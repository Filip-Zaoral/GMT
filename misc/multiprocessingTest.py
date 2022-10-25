
# Definitions:
p = 8
N = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
a = 5
b = 10

# Function to execute in parallel:
def f(a, b, n):
    from time import sleep
    sleep(5)
    r = n + a * b
    s = ["Screw you GIL, KEKW."]
    return r,s

# We execute the function in parallel using Pool.map()
def g():
    from multiprocessing import Pool
    from functools import partial
    from time import perf_counter
    T = 0
    for i in range(1):
        t = perf_counter()
        with Pool(p) as P:
            result,L = zip(*P.map(partial(f, a, b), N))
            print(list(result))
            print(list(L))
            print(P._processes)
        t = perf_counter() - t
        T +=  t
    T = T/1
    print(T)
