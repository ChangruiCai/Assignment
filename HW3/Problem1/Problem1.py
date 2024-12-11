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


def loss(y, y_pred):
    return np.mean((y - y_pred) ** 2)


def grad_m(y, y_pred, x):
    return -2 * np.mean(x * (y - y_pred))


def grad_b(y, y_pred):
    return -2 * np.mean(y - y_pred)



def linear(epochs, lr):
    m = np.random.rand()
    b = np.random.rand()

    for epoch in range(epochs):
        y_pred = m * x_real + b

        loss_value = loss(y_real, y_pred)

        dm = grad_m(y_real, y_pred, x_real)
        db = grad_b(y_real, y_pred)

        m -= lr * dm
        b -= lr * db

        if epoch % 100 == 0:
            print(f'Epoch:{epoch}, Loss:{loss_value}')

    y_pred = m * x_real + b
    plt.figure()
    plt.plot(x_real, y_real, alpha=0.5, label='Actual Data(Noisy)')
    plt.plot(x_real, y_pred, 'r--', label='Predicted Data(Linear)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Linear Model')
    plt.legend()
    plt.show()


linear(10000, 0.001)
