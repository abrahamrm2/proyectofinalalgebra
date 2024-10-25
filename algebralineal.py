import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt

def calcular_inversa():
    try:
        matriz = np.array(eval(entry_matriz.get()))
        if matriz.shape[0] != matriz.shape[1]:
            raise ValueError("La matriz debe ser cuadrada.")
        inversa = np.linalg.inv(matriz)
        messagebox.showinfo("Inversa", f"La inversa es:\n{inversa}")
    except np.linalg.LinAlgError:
        messagebox.showerror("Error", "La matriz no es invertible.")
    except ValueError as ve:
        messagebox.showerror("Error", str(ve))
    except Exception as e:
        messagebox.showerror("Error", f"Entrada no válida: {str(e)}")

def multiplicar_matrices():
    try:
        matriz1 = np.array(eval(entry_matriz1.get()))
        matriz2 = np.array(eval(entry_matriz2.get()))
        resultado = np.dot(matriz1, matriz2)
        messagebox.showinfo("Resultado", f"El resultado es:\n{resultado}")
    except ValueError as ve:
        messagebox.showerror("Error", str(ve))
    except Exception as e:
        messagebox.showerror("Error", f"Entrada no válida: {str(e)}")

def resolver_ecuaciones():
    try:
        matriz_a = np.array(eval(entry_matriz_a.get()))
        matriz_b = np.array(eval(entry_matriz_b.get()))
        if matriz_a.shape[0] != matriz_b.shape[0] or matriz_a.shape[1] != matriz_b.shape[0]:
            raise ValueError("Las dimensiones de las matrices no coinciden.")

        if metodo_var.get() == "gauss_jordan":
            solucion = gauss_jordan(matriz_a, matriz_b)
        else:
            solucion = cramer(matriz_a, matriz_b)

        # Convertir la solución a una lista de números flotantes
        solucion_legible = [float(sol) for sol in solucion]
        messagebox.showinfo("Solución", f"La solución es:\n{solucion_legible}")
        graficar_ecuaciones(matriz_a, matriz_b)
    except np.linalg.LinAlgError:
        messagebox.showerror("Error", "El sistema no tiene solución o tiene infinitas soluciones.")
    except ValueError as ve:
        messagebox.showerror("Error", str(ve))
    except Exception as e:
        messagebox.showerror("Error", f"Entrada no válida: {str(e)}")


def gauss_jordan(matriz_a, matriz_b):
    # Implementación del método de Gauss-Jordan
    n = len(matriz_b)
    augmented_matrix = np.hstack((matriz_a, matriz_b.reshape(-1, 1)))
    for i in range(n):
        augmented_matrix[i] = augmented_matrix[i] / augmented_matrix[i][i]
        for j in range(n):
            if i != j:
                augmented_matrix[j] -= augmented_matrix[i] * augmented_matrix[j][i]
    return augmented_matrix[:, -1]

def cramer(matriz_a, matriz_b):
    # Implementación del método de Cramer
    det_a = np.linalg.det(matriz_a)
    if det_a == 0:
        raise np.linalg.LinAlgError("El determinante es cero, el sistema no tiene solución única.")
    
    n = matriz_a.shape[0]
    soluciones = []
    for i in range(n):
        matriz_temp = matriz_a.copy()
        matriz_temp[:, i] = matriz_b
        soluciones.append(np.linalg.det(matriz_temp) / det_a)
    return soluciones

def graficar_ecuaciones(matriz_a, matriz_b):
    # Suponiendo que la matriz A tiene 2 ecuaciones y 2 incógnitas
    x = np.linspace(-10, 10, 400)
    y1 = (matriz_b[0] - matriz_a[0, 0] * x) / matriz_a[0, 1]
    y2 = (matriz_b[1] - matriz_a[1, 0] * x) / matriz_a[1, 1]

    plt.figure()
    plt.plot(x, y1, label='Ecuación 1')
    plt.plot(x, y2, label='Ecuación 2')
    plt.axhline(0, color='black', lw=0.5, ls='--')
    plt.axvline(0, color='black', lw=0.5, ls='--')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.title('Gráfico de Ecuaciones')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.legend()
    plt.show()

# Crear la ventana principal
root = tk.Tk()
root.title("Calculadora de Álgebra Lineal")

# Entradas para la inversa de la matriz
tk.Label(root, text="Ingrese la matriz para calcular la inversa:").pack()
entry_matriz = tk.Entry(root)
entry_matriz.pack()
tk.Button(root, text="Calcular Inversa", command=calcular_inversa).pack()

# Entradas para multiplicación de matrices
tk.Label(root, text="Ingrese la primera matriz:").pack()
entry_matriz1 = tk.Entry(root)
entry_matriz1.pack()
tk.Label(root, text="Ingrese la segunda matriz:").pack()
entry_matriz2 = tk.Entry(root)
entry_matriz2.pack()
tk.Button(root, text="Multiplicar Matrices", command=multiplicar_matrices).pack()

# Entradas para resolver sistemas de ecuaciones
tk.Label(root, text="Ingrese la matriz A (coeficientes):").pack()
entry_matriz_a = tk.Entry(root)
entry_matriz_a.pack()
tk.Label(root, text="Ingrese la matriz B (resultados):").pack()
entry_matriz_b = tk.Entry(root)
entry_matriz_b.pack()

# Opción para seleccionar el método
metodo_var = tk.StringVar(value="gauss_jordan")
tk.Radiobutton(root, text="Método de Gauss-Jordan", variable=metodo_var, value="gauss_jordan").pack()
tk.Radiobutton(root, text="Método de Cramer", variable=metodo_var, value="cramer").pack()

tk.Button(root, text="Resolver Ecuaciones", command=resolver_ecuaciones).pack()

# Ejecutar la aplicación
root.mainloop()
