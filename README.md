# Алгоритм Витерби мягкого декодирования линейных блоковых кодов

Имя входного файла: input.txt
Имя выходного файла: output.txt
Ограничение по времени: 15 секунд
Ограничение по памяти: 256 мегабайт
Необходимо реализовать кодер линейного блокового кода и его декодер на основе алгоритма Ви-
терби. Исходными данными являются длинаn, размерностьkи порождающая матрицаGдвоичного
линейного блокового кода.

## Форматвходныхданных

файл input.txt должен начинаться следующим образом:
n k
G
Далее должны быть представлены команды.

## Форматвыходныхданных

Первой строкой в файле output.txt должна быть последовательность|Vi|, 06 i 6 n, где|Vi|
число узлов в решетке на ярусеi. Далее должны быть приведены результаты выполнения команд.
Моделирование следует производить для случая канала с двоичной амплитудно-импульсной моду-
ляцией и аддитивным белым гауссовским шумом. Под уровнем шума следует понимать отношение
сигнал/шум на бит, выраженное в децибелах.

## Пример

```
input.txt output.txt
8 4
1 1 1 1 1 1 1 1
1 1 1 1 0 0 0 0
1 1 0 0 1 1 0 0
1 0 1 0 1 0 1 0
Encode 1 0 0 0
Decode -1.0 1.0 1 1 1 1 1 1.
Simulate 3 100000 100
Simulate 4 100000 100
```
### 1 2 4 8 4 8 4 2 1

### 1 1 1 1 1 1 1 1

### 0 0 0 0 0 0 0 0

### 2.56E-

### 9.31E-
