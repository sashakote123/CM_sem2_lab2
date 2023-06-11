import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Ваши данные точек
X = [1, 2, 3, 4, 5]
Y = [2, 4, 6, 8, 10]
Z = [3, 6, 9, 12, 15]

# Преобразование данных точек в двумерные массивы
points = np.array([X, Y, Z])

# Аппроксимация плоскости
coefficients = np.polyfit(points[0], points[1], 1)
a, b = coefficients

# Создание сетки точек для плоскости
x_range = np.linspace(min(X), max(X), 10)
y_range = np.linspace(min(Y), max(Y), 10)
X_grid, Y_grid = np.meshgrid(x_range, y_range)
Z_plane = a * X_grid + b * Y_grid

# Создание трехмерного графика
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Построение точек
ax.scatter(X, Y, Z, c='r', marker='o')

# Построение плоскости
ax.plot_surface(X_grid, Y_grid, Z_plane, alpha=0.5)

# Настройка осей
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Отображение графика
plt.show()
