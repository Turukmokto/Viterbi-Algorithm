#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include <map>
#include <random>

#define _CRT_SECURE_NO_WARNINGS

#pragma warning (disable : 4996)

using namespace std;

int n, k;
random_device rd;
mt19937 gen(rd());
vector<bool> used_rows;
vector<double> sizes_of_tiers;
vector<vector<int>> G, G_span, count_vertex;
vector<pair<int, int>> active_segments;
struct node;
vector<vector<node>> graph;

bool check_line(int i, int j);

bool is_a_rib(int cur_tier, int v, int u);

double gen_normal(double stddev);

double simulate(int noise, int count_tests, int count_mis);

void in(vector<double> &arr, int sz);

void out(vector<double> &arr);

void line_addition(int cur_column, int cur_line);

void swap_lines(int last_active_line, int cur_line);

void forward_transform();

void backward_transform();

void matrix_to_span();

void make_count_vertex();

void addition_mask_values(int cur_tier, int cur_ver, const vector<int> &lines_values);

void make_mask(int mask, vector<int> &lines_values);

void make_inf_bits();

void add_inform(int cur_tier, int v, set<pair<int, int>> &un);

void paint_rib(int cur_tier, int v, int u, int color);

void make_ribs();

void make_graph();

void clean();

void relax(const vector<double> &code, int cur_tier, int v, int type, int next);

vector<double> path_restoration();

vector<double> decode(vector<double> &code);

vector<double> noise_row(vector<double> &r, double noiseLevel);

vector<double> stdRow();

vector<double> encode(vector<double> &word);


vector<double> encode(vector<double> &word) { //Кодировка слова
    vector<double> encode_ans(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            // Умножение вектора на порождающую матрицу
            encode_ans[i] = (int) (word[j] * G[j][i] + encode_ans[i]) % 2;
        }
    }
    return encode_ans;
}

bool check_line(int i, int j) { // Проверка активированности строки в матрице
    return !used_rows[j] && G_span[j][i] == 1;
}


void line_addition(int cur_column,int cur_line) { // Поиск строк в матрице, которые одновременно активириуются или деактивируются, и их устранение
    for (int j = 0; j < k; j++) {
        if (check_line(cur_column,j)) { // Если на соответствующем столбце есть еще одна активированная строка - прибавим к ней стартовую
            for (int i = 0; i < n; i++) {
                G_span[j][i] = (G_span[j][i] + G_span[cur_line][i]) % 2;
            }
        }
    }
}

void swap_lines(int last_active_line, int cur_line) { // Меняем строки и информацию о них местами
    swap(G_span[cur_line], G_span[last_active_line]);
    swap(used_rows[cur_line], used_rows[last_active_line]);
    swap(active_segments[cur_line], active_segments[last_active_line]);
}


void
forward_transform() { // Первый этап превращения матрицы в минимальную спэнову форму, сделаем так, чтобы начала активации строк попарно не совпадали
    int last_active_line = 0; // Место, где должна стоять последняя активированная строка, чтобы матрица имела диагональный вид
    used_rows.assign(k, false);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            if (check_line(i, j)) { // Находим первую строчку, которая активируется в этом столбце
                used_rows[j] = true;
                active_segments[j].first = i; // Запоминаем момент ативации
                line_addition(i, j); // Добавляем строку к тем, которые также активировались в этом столбце
                if (j != last_active_line) {
                    swap_lines(last_active_line, j); // Если строка не стоит на нужном месте - меняем ее местами
                }
                last_active_line++;
                break;
            }
        }
    }
}


void
backward_transform() { // Второй этап превращения матрицы в минимальную спэнову форму, сделаем так, чтобы концы активации строк попарно  не совпадали
    used_rows.assign(k, false);
    for (int i = n - 1; i >= -1; i--) {
        for (int j = k - 1; j > -1; j--) {
            if (check_line(i, j)) { // Находим первую строчку, которая деактивируется в этом столбце
                used_rows[j] = true;
                active_segments[j].second = i; // Запоминаем момент деактивации
                line_addition(i, j); // Добавляем строку к тем, которые также деактивировались в этом столбце
                break; // Так как матрица уже диагональная, то после обновлений мы не изменим начала активации строк
            }
        }
    }
}


void matrix_to_span() { // Превращение матрицы в минимальную спэнову форму
    forward_transform();
    backward_transform();
}


void make_count_vertex() { // Запомним номера столбцов относящихся к конкретным слоям графа
    for (int i = 0; i < k; i++) {
        for (int j = active_segments[i].first + 1; j <= active_segments[i].second; j++) {
            count_vertex[j].push_back(i);
            sizes_of_tiers[j] *= 2; // А также посчитаем будущее количество вершин в каждом слое
        }
    }
}

struct node { // Вершина графа
    int rib_0 = -1, rib_1 = -1; // Индексы перехода по ребрам 0 и 1
    vector<pair<int, int>> inf_bits; // Информация об активных столбцах, соответствующих вершине, а так же об их значениях
    pair<int, int> parent = {-1,-1}; // Вершина на предыдущем шаге, из которой был оптимальный путь, а так же цвет последнего ребра
    double dist = 1e9; // Минимальная разница в сообщении, с которой можно добраться до вершины
};

void addition_mask_values(int cur_tier, int cur_ver,const vector<int> &lines_values) { // Заполнение информации об активированных столбцах и их значениях в вершинах графа
    for (int i = 0; i < count_vertex[cur_tier].size(); i++) {
        graph[cur_tier][cur_ver].inf_bits.emplace_back(count_vertex[cur_tier][i], lines_values[i]);
    }
}

void make_mask(int mask, vector<int> &lines_values) { // Генерация битовых значений для вершины посредством маски
    int ind = 0;
    while (mask > 0) {
        if (mask % 2 == 1) {
            lines_values[ind] = 1;
        }
        ind++;
        mask /= 2;
    }
}

void
make_inf_bits() { // Создание слоев и вершин в графе и их инициализация в зависимости от соответствующих слоям столбцов
    for (int i = 0; i < n + 1; i++) {
        graph[i].resize((int) sizes_of_tiers[i]);
        for (int j = 0; j < sizes_of_tiers[i]; j++) {
            vector<int> columns_values(count_vertex[i].size(), 0); // Значения столбца в конкретной вершине слоя
            make_mask(j, columns_values); // Генерируем все возможные варианты значений для конкретного слоя
            addition_mask_values(i, j, columns_values); // Обновляем информацию
        }
    }
}


bool is_a_rib(int cur_tier, int v, int u) { // Проверка существования ребра между двумя вершинами двух соседних слоев
    for (int l1 = 0; l1 < graph[cur_tier][v].inf_bits.size(); l1++) {
        for (int l2 = 0; l2 < graph[cur_tier + 1][u].inf_bits.size(); l2++) {
            // Если у вершин имеется один и тот же активированный столбец и он имеет в них разные значения => ребра между ними нет
            if (graph[cur_tier][v].inf_bits[l1].first == graph[cur_tier + 1][u].inf_bits[l2].first &&
                graph[cur_tier][v].inf_bits[l1].second != graph[cur_tier + 1][u].inf_bits[l2].second) {
                return false;
            }
        }
    }
    return true;
}


void add_inform(int cur_tier, int v, set<pair<int, int>> &un) { // Обьединение информации по столбцам и их значениям
    for (auto &inf_bit : graph[cur_tier][v].inf_bits) {
        un.insert(inf_bit);
    }
}

void paint_rib(int cur_tier, int v, int u, int color) { // Проведение нужного ребра в зависимости от цвета
    if (color == 1) {
        graph[cur_tier][v].rib_1 = u;
    } else {
        graph[cur_tier][v].rib_0 = u;
    }
}

void make_ribs() { // Создание ребер в графе и их дальнейшая покраска
    for (int cur_tier = 0; cur_tier < n; cur_tier++) {
        for (int v = 0; v < graph[cur_tier].size(); v++) {
            for (int u = 0; u < graph[cur_tier + 1].size(); u++) {
                bool check = is_a_rib(cur_tier, v, u); // Проверка на существование ребра
                if (check) {
                    set<pair<int, int>> un; // Выпишем в вектор все значения, соответствующие столбцам рассматриваемых вершин
                    add_inform(cur_tier, v, un);
                    add_inform(cur_tier + 1, u, un);
                    int color = 0;
                    for (auto it : un) { // Чтобы узнать цвет ребра, перемножим значения активных столбцов в вершинах и спеновой матрице
                        color += it.second * G_span[it.first][cur_tier];
                        color %= 2;
                    }
                    paint_rib(cur_tier, v, u, color); // Покрасим соответствующее ребро
                }
            }
        }
    }
}


void make_graph() { // Создание графа
    make_inf_bits();
    make_ribs();
}


void clean() { // Подготовка графа к новому запуску
    for (int cur_tier = 0; cur_tier < n + 1; cur_tier++) {
        for (auto &v : graph[cur_tier]) {
            v.dist = 1e9;
            v.parent = {-1, -1};
        }
    }
    graph[0][0].dist = 0;
}

void relax(const vector<double> &code, int cur_tier, int v, int type, int next) { // Обновление значения расстояния в вершине, если нашли более оптимальный путь
    if (next != -1) {
        if (graph[cur_tier + 1][next].dist > graph[cur_tier][v].dist + abs(code[cur_tier] + (-1.0 + 2 * type))) {
            graph[cur_tier + 1][next].dist = graph[cur_tier][v].dist + abs(code[cur_tier] + (-1.0 + 2 * type));
            graph[cur_tier + 1][next].parent = {v, type};
        }
    }
}

vector<double> path_restoration() { // Восстановление оптимального пути
    int cur_v = 0, cur_sl = n;
    vector<double> deans;
    while (cur_sl != 0) {
        deans.push_back(
                graph[cur_sl][cur_v].parent.second); // Так как мы запоминали не только предка но и цвет последнего ребра - восстановить ответ много легче
        cur_v = graph[cur_sl][cur_v].parent.first;
        cur_sl--;
    }
    reverse(deans.begin(), deans.end()); // Путь восстановлен в обратном порядке, поэтому развернем его
    return deans;
}

vector<double> decode(vector<double> &code) { // Расшифровка заданного слова
    clean();
    for (int cur_tier = 0; cur_tier < n; cur_tier++) {
        for (int v = 0; v < graph[cur_tier].size(); v++) {
            relax(code, cur_tier, v, 0, graph[cur_tier][v].rib_0); // Пытаемся улучшить путь используя наши 2 ребра
            relax(code, cur_tier, v, 1, graph[cur_tier][v].rib_1);
        }
    }
    return path_restoration(); // Восстанавливаем и выводим путь, то есть расшифрованное сообщение
}


double gen_normal(double stddev) { // Генерация нормальных помех
    normal_distribution<> X(0, stddev);
    return X(gen);
}

vector<double> noise_row(vector<double> &r, double noiseLevel) { // Генерация закодированного сообщения с помехами
    vector<double> noised(r.size());

    double stddev = sqrt(0.5 * pow(10, -noiseLevel / 10) * (double) n / (double) k);
    for (int i = 0; i < r.size(); ++i) {
        noised[i] = 1.0 - 2.0 * r[i] + gen_normal(stddev);
    }

    return noised;
}

vector<double> stdRow() { // Генерация обычного сообщения без помех
    vector<double> r(k);
    for (int i = 0; i < k; i++) {
        r[i] = gen_normal(1) < 0;
    }
    return r;
}


double simulate(int noise, int count_tests, int count_mis) { // Имитация принятия множества сообщений и их декодирование
    int cur_mis = 0; // Сколько тестов провалились
    int cur_test; // Текущий тест
    for (cur_test = 0; cur_test < count_tests; cur_test++) {
        // Генерация обычного сообщения, сообщения с помехами и их декодирование
        vector<double> std_row = stdRow();
        vector<double> std_code = encode(std_row);
        vector<double> noise_code = noise_row(std_code, noise);
        vector<double> deans = decode(noise_code);
        if (deans != std_code) {
            cur_mis++;
        }
        if (cur_mis == count_mis) { // Если количество допутимых ошибок превышено - прекращаем работу
            break;
        }
    }
    cur_test++;
    return (double) cur_mis / (double) cur_test; // Выводим процент ошибок относительно общего числа тестов
}

void in(vector<double> &arr, int sz) { // Ввод массива
    for (int i = 0; i < sz; i++) {
        cin >> arr[i];
    }
}

void out(vector<double> &arr) { // Вывод массива
    for (auto it : arr) {
        cout << it << " ";
    }
    cout << endl;
}


int main() {
    ios_base::sync_with_stdio(false);
    freopen("input.txt", "r", stdin);
    //freopen ("output.txt", "w", stdout);
    cin.tie(nullptr);
    cout.tie(nullptr);
    cin >> n >> k;
    //Задаем нужные стартовые размеры обьектов
    active_segments.resize(k, {-1, -1});
    G.resize(k, vector<int>(n));
    G_span.resize(k, vector<int>(n));
    count_vertex.resize(n + 1);
    sizes_of_tiers.resize(n + 1, 1);
    graph.resize(n + 1);
    //Считываем порождающую матрицу
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            cin >> G[i][j];
            G_span[i][j] = G[i][j];
        }
    }
    matrix_to_span();
    make_count_vertex();
    out(sizes_of_tiers);
    make_graph();
    string message;
    // Принимаем запросы
    while (cin >> message) {
        if (message == "Encode") {
            vector<double> word(k), encode_ans;
            in(word, k);
            encode_ans = encode(word);
            out(encode_ans);
        }
        if (message == "Decode") {
            vector<double> code(n), deans;
            in(code, n);
            deans = decode(code);
            out(deans);
        }
        if (message == "Simulate") {
            int noise, count_tests, count_mis;
            cin >> noise >> count_tests >> count_mis;
            cout << simulate(noise, count_tests, count_mis) << endl;
        }
    }
    return 0;
}