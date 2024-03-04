В директории singular point решение задачи Коши для системы дифференциальных уравнений явным методом Рунге-Кутты 4-го порядка. Для разработки интерфейса программы использовался Qt. Интерфейс представляет собой открывающееся главное окно, в котором пользователь вводит данные, а затем выбирает одну из трех опций для отрисовки графика (y1(t), y2(t), y2(y1)) и по нажатию кнопки «нарисовать» открывается в зависимости от выбора опции другое окно  и закрывается предыдущее, если была выбрана первая или вторая опция, то в новом окне показывается соответствующий график функции, далее можно нажать на кнопки либо для выхода из программы либо для возвращения к главному окну для ввода новых данных, если была выбрана третья опция, то также отображается соответствующий график функции, а также название особой точки и координаты стационарной точки, в этом окне также аналогично можно либо выйти из программы либо вернуться к главному окну. Для запуска исполняемого файла task в системе должен стоять gnuplot, так как с помощью него отрисовываются графики.
В директории solvers различные явные, неявные методы решения задачи Коши для системы дифференциальных уравнений.
