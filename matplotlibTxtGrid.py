import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

# Чтение узлов из файла
nodes = []
with open('nodes.txt', 'r') as file:
    for line in file:
        x, y = map(float, line.strip().split(','))
        nodes.append([x, y])

# Преобразуем узлы в numpy-массив
nodes = np.array(nodes)

# Выполняем триангуляцию
tri = Delaunay(nodes)

# Устанавливаем пороги для фильтрации
min_angle = 20  # Минимальный угол в градусах
max_angle = 120  # Максимальный угол в градусах


# Функция для вычисления углов треугольника
def calculate_angles(triangle):
    p1, p2, p3 = nodes[triangle]
    # Векторы сторон треугольника
    a = p2 - p1
    b = p3 - p2
    c = p1 - p3

    # Длины сторон
    len_a = np.linalg.norm(a)
    len_b = np.linalg.norm(b)
    len_c = np.linalg.norm(c)

    angle1 = np.degrees(np.arccos(np.clip(np.dot(-a, c) / (len_a * len_c), -1.0, 1.0)))
    angle2 = np.degrees(np.arccos(np.clip(np.dot(-b, a) / (len_b * len_a), -1.0, 1.0)))
    angle3 = 180 - angle1 - angle2
    return angle1, angle2, angle3


filtered_triangles_by_angle = []
for triangle in tri.simplices:
    angles = calculate_angles(triangle)
    if all(min_angle <= angle <= max_angle for angle in angles):
        filtered_triangles_by_angle.append(triangle)

# Визуализация результатов
plt.figure(figsize=(12, 6))


plt.subplot(1, 2, 1)
plt.triplot(nodes[:, 0], nodes[:, 1], filtered_triangles_by_angle, color='blue')
plt.plot(nodes[:, 0], nodes[:, 1], 'o', color='red', markersize=0.5)
plt.axis('equal')
plt.grid(False)

plt.show()

# Сохранение отфильтрованных треугольников в файлы
with open('elements.txt', 'w') as file:
    for triangle in filtered_triangles_by_angle:
        file.write(f"{triangle[0]}, {triangle[1]}, {triangle[2]}\n")

