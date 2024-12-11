import numpy as np
import matplotlib.pyplot as plt

# Seed for reproducibility
np.random.seed(42)

# Generate 10 random x values within a range
N = 10
x_real = np.linspace(0, 5, N)

# Parameters for the function (can use the previously fitted values or set randomly)
n_true = 0.06
a_true = 0.25
m_true = 0.57
b_true = 0.11

# Generate corresponding y values based on the function with added noise
noise = 0.001 * np.random.normal(0, 0.1, size=x_real.shape)

# Add Gaussian noise
y_real = n_true * np.exp(-a_true * (m_true * x_real + b_true) ** 2) + noise
    
def grad_n(y, y_pred, y_int, a):
    return -2 * np.mean((y - y_pred) * np.exp(-a * y_int))

def grad_a(y, y_pred, y_int):
    return 2 * np.mean((y - y_pred) * y_pred * (-y_int))

def grad_m(y, y_pred, y_base, x, a):
    return 2 * np.mean((y - y_pred) * y_pred * (-a) * 2 * y_base * x)

def grad_b(y, y_pred, y_base, a):
    return 2 * np.mean((y - y_pred) * y_pred * (-a) * 2 * y_base)

def nonlinear(epochs, lr):
    #seed
    n = (np.random.rand() + np.random.rand())/2
    a = np.random.rand()
    m = np.random.rand()
    b = np.random.rand()

    for epoch in range(epochs):
        y_base = m * x_real + b
        y_int = y_base ** 2
        #mse
        y_pred = n * np.exp(-a * y_int)

        

        dn = grad_n(y_real, y_pred, y_int, a)
        da = grad_a(y_real, y_pred, y_int)
        dm = grad_m(y_real, y_pred, y_base, x_real, a)
        db = grad_b(y_real, y_pred, y_base, a)

        n -= lr * dn
        a -= lr * da
        m -= lr * dm
        b -= lr * db



    y_int = (m * x_real + b) ** 2
    y_pred = n * np.exp(-a * y_int)

    plt.figure()
    plt.plot(x_real, y_real, label='Actual Data')
    plt.plot(x_real, y_pred, 'r--', label='Predicted Data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Nonlinear Model')
    plt.legend()
    plt.show()

#print(f'Epoch:{epoch}, Loss:{loss_value}')
nonlinear(10000, 0.001)

