Программа реализует Метод Конечных Элементов для решения дифференциальных уравнений классической теории упругости.
Решается задача Кирша, для четверти пластины с отверстием.
Сетка строится в модуле gridGeneration.cpp, при построении сетки используется библиотека gmsh. Для ручного контроля построения сетки используется модуль matplotlibTxtGrid.py в котором строится триангуляция Делоне, и отрисовка получившейся конечно-элементной треугольной сетки. 
После проведения расчетов в основном модуле MKEvar0.cpp данные передаются в Python модуль ResultsGraph.py для проверки корректности решения - отображения графиков сравнения численного и аналитического решения.
