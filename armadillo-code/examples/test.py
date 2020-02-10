import numpy as np




def numerical_extrapolation(H_0, H_t, N=100):
    def H(n):
        return H_0 + ((H_t - H_0)/N * n)

    R = 1.0
    for i in range(N):
        R *= pow(N, (1/N)) * np.exp(-0.5*(H(i) + H(i+1)))
    return R

def analytical_extrapolation(H_0, H_t, N=100):
    return  N*np.exp(-0.5*(H_0 + H_t) - ((N-1)/2) * (H_t + H_0))


for _ in range(100):
    H_0 = np.random.random()
    H_t = np.random.random()
    print(numerical_extrapolation(H_0, H_t))
    print(analytical_extrapolation(H_0, H_t))
    assert np.round(numerical_extrapolation(H_0, H_t), 10) == np.round(analytical_extrapolation(H_0, H_t), 10)
