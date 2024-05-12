//=======================================================================
// Задача коммивояжера (TSP) точное решение — метод ветвей и границ (by ReBuilder)
// https://habr.com/ru/post/708072/
//=======================================================================
const unsigned int inf = 0x7fffffff;  // 2147483647

typedef struct TVector {
    unsigned int src;
    unsigned int dest;
} TVector;

typedef struct TMatrixObj {
    unsigned int global_min;  // Глобальный текущий минимум расчёта
    unsigned int n;           // Реальный размер матрицы
    unsigned int step;        // Шаг алгоритма
    TVector *v;               // Указатель на массив векторов сохранённых рёбер на текущем шаге [откуда, куда]
    TVector *v_best;          // Указатель на массив векторов, лучший на текущий момент [откуда, куда]
} TMatrixObj;

void matrix_evaluation(struct TMatrixObj *m, unsigned int *pm, unsigned int index, unsigned int n);
//=======================================================================
// Сохраняет два минимальных элемента при добавлении нового
void two_lows(unsigned int *a, unsigned int *b, unsigned int c)
{
    if (c <= *a)
    {
        *b = *a;
        *a = c;
    }
    else
        if (c < *b)
            *b = c;
}
//=======================================================================
void matrix_reduction(unsigned int *pm, unsigned int n, unsigned int *max_sum, unsigned int *max_i, unsigned int *max_j)
{
    unsigned int sum = 0;

    // Организуем кеш для подсчёта веса строки
    unsigned int line_cache[n + 1];

    // Для каждой строки находим минимум
    for (unsigned int i = 1; i < n; i++)
    {
        unsigned int min1 = pm[i * n + 1];
        unsigned int min2 = inf;
        two_lows(&min1, &min2, pm[i * n + 2]);

        for (unsigned int j = 3; j < n; j++)
            two_lows(&min1, &min2, pm[i * n + j]);

        // Если хоть одна точка источник не имеет выходов то обрываем расчёт
        if (min1 == inf)
        {
            pm[0] = inf;
            return;
        }

        // Вычитаем минимум из каждой строки (Редукция строк)
        if (min1 > 0)
        {
            if (min2 < inf)
                min2 -= min1;

            for (unsigned int j = 1; j < n; j++)
            if (pm[i * n + j] < inf)
                pm[i * n + j] -= min1;
            // Прибавляем минимальный элемент к нижней границе
            sum += min1;
        }
        line_cache[i] = min2;
    }

    *max_sum = 0;
    *max_i = inf;

    // Для каждого столбца находим минимум
    for (unsigned int i = 1; i < n; i++)
    {
        unsigned int min1 = pm[1 * n + i];
        unsigned int min2 = inf;
        two_lows(&min1, &min2, pm[2 * n + i]);

        for (unsigned int j = 3; j < n; j++)
            two_lows(&min1, &min2, pm[j * n + i]);

        // Если хоть одна точка назначения не имеет входов то обрываем расчёт
        if (min1 == inf)
        {
            pm[0] = inf;
            return;
        }

        if (min1 > 0)
        {
            if (min2 < inf)
                min2 -= min1;

            // Вычитаем минимум из каждого столбца (Редукция столбцов)
            // Находим элемент для разбиения, и верхнюю оценку
            unsigned int temp_val;
            for (unsigned int j = 1; j < n; j++)
            {
                temp_val = pm[j * n + i];
                if (temp_val < inf)
                {
                    temp_val -= min1;
                    pm[j * n + i] = temp_val;
                    if (temp_val == 0)
                    {
                        temp_val = line_cache[j] + min2;
                        if (temp_val >= *max_sum)
                        {
                            *max_sum = temp_val;
                            *max_i = j;
                            *max_j = i;
                        }
                    }
                }
            }

            // Прибавляем минимальный элемент к нижней границе
            sum += min1;

        }
        else
        {
            // Находим элемент для разбиения, и верхнюю оценку
            unsigned int temp_val;
            for (unsigned int j = 1; j < n; j++)
                if (pm[j * n + i] == 0)
                {
                    temp_val = line_cache[j] + min2;
                    // выбираем значение с максимальной суммой минимумов по строке и столбцу
                    if (temp_val >= *max_sum)
                    {
                        *max_sum = temp_val;
                        *max_i = j;
                        *max_j = i;
                    }
                }
        }
    }

    // Нижняя граница - стоимость меньше которой невозможно построить маршрут
    // Если нет элемента для разбиения, то блочим ветку
    if (*max_i == inf)
        pm[0] = inf;
    else
        pm[0] += sum;
}
//=======================================================================
int head_search(unsigned int *pm, struct TVector *v, unsigned int n, unsigned int index, unsigned int dj)
{
    unsigned int l = 0;
    do
    {
        if (v[l].src == dj)
        {
            dj = v[l].dest;
            l = 0;
            continue;
        }
        l++;
    } while (l < index);

    for (l = 1; l < n; l++)
        if (pm[l * n + 0] == dj)
            break;
    return l;
}
//=======================================================================
int tail_search(unsigned int *pm, struct TVector *v, unsigned int n, unsigned int index, unsigned int di)
{
    unsigned int l = 0;
    do
    {
        if (v[l].dest == di)
        {
            di = v[l].src;
            l = 0;
            continue;
        }
        l++;
    } while (l < index);

    for (l = 1; l < n; l++)
        if (pm[0 * n + l] == di)
            break;
    return l;
}
//=======================================================================
void martix_in_depth(struct TMatrixObj *m, unsigned int *pm, unsigned int n, unsigned int index, unsigned int di, unsigned int dj)
{
    unsigned int new_matrix[n * n];
    // Копируем матрицу в массив меньшего размера
    {
        unsigned int *src = pm;
        unsigned int *dest = new_matrix;
        for (unsigned int i = 0; i <= n; i++)
        {
            if (i == di)
            {
                src += n + 1;
                continue;
            }
            for (unsigned int j = 0; j <= n; j++)
            {
                if (j == dj)
                {
                    src++;
                    continue;
                }
                *dest++ = *src++;
            }
        }
    }

    if (n > 3)
    {
        // Запрещаем переходы в уже пройденные узлы графа.
        // Чтобы избежать подциклы в векторе найденых рёбер ищем голову и хвост цепочек,
        // вместе с которыми добавленное ребро образует единую ветвь.
        // Для найденый индексов в матрице меньшего размера ищем обратные индексы
        unsigned int i = head_search(new_matrix, m->v, n, index, m->v[index].dest);
        unsigned int j = tail_search(new_matrix, m->v, n, index, m->v[index].src);
        // Зануляем обратный элемент в матрице меньшего размера.
        // Должен находиться всегда!!!
        new_matrix[i * n + j] = inf;
    }
    else if (n == 2)
    {
        // Если в матрице остался только один элемент
        new_matrix[0] += new_matrix[3];

        // Если текущее решение лучшее из найденных, запоминаем его
        if (new_matrix[0] < m->global_min)
        {
            m->global_min = new_matrix[0];
            for (unsigned int i = 0; i <= index; i++)
            {
               m->v_best[i].src = m->v[i].src;
               m->v_best[i].dest = m->v[i].dest;
            }
            //memcpy(m->v_best, m->v, (m->n-1) * sizeof(TVector));
            m->v_best[index + 1].src = new_matrix[1 * n + 0];
            m->v_best[index + 1].dest = new_matrix[0 * n + 1];
        }
        return;
    }
    // Производим раcчёт для слоя ниже
    matrix_evaluation(m, new_matrix, index + 1, n);
}
//=======================================================================
void matrix_evaluation(struct TMatrixObj *m, unsigned int *pm, unsigned int index, unsigned int n)
{
    unsigned int max_sum, max_i, max_j;
    do
    {
        m->step++;
        // Производим приведение матрицы, находим элемент для разбиения, и верхнюю оценку
        matrix_reduction(pm, n, &max_sum, &max_i, &max_j);
        // Если нет возможного пути, то блочим ветку
        if (pm[0] >= inf)
            return;

        // Расчитываем множество М1
        // Если лучший из найденых пока маршрутов больше нижней границы,
        // то вычёркиваем текущую строку и столбец и переходим к уменьшиной матрице
        if (m->global_min > pm[0])
        {
            // Сохраняем индексы элемента разбиения в вектор локального результата
            m->v[index].src = pm[max_i * n + 0];
            m->v[index].dest = pm[0 * n + max_j];
            // Погружаемся на слой ниже
            martix_in_depth(m, pm, n - 1, index, max_i, max_j);
        }
        // Расчитываем множество М2
        // Исключаем пункт из множества
        pm[max_i * n + max_j] = inf;
        // Производим повторную оценку после вычерка вершины (множество М2)
        // Так как хвостовая рекурсия дороже итерации, оформляем цикл
    }
    while (m->global_min >= pm[0] + max_sum);
}
//=======================================================================


#ifdef __cplusplus
extern "C" {
#endif
    int tsp_branch(unsigned int n, int *p)
    {
//        if ((n < 2) || (n > 64))
//           return 0;

        struct TMatrixObj m;
        m.n = n;
        struct TVector v_best[n];
        struct TVector v[n - 1];
        m.v = v;
        m.v_best = v_best;
        unsigned int work_matrix[(n + 1) * (n + 1)];

        //-----------------------------------------------------------------------
        // Производим первичную инициализацию рабочей матрицы
        m.global_min = inf;
        m.step = 0;
        work_matrix[0] = 0;
        for (unsigned int i = 1; i < m.n + 1; i++)
        {
            work_matrix[i * (m.n + 1) + 0] = i;
            work_matrix[0 * (m.n + 1) + i] = i;
        }
        {
            int *src = p;
            for (unsigned int i = 1; i < m.n + 1; i++)
                for (unsigned int j = 1; j < m.n + 1; j++, src++)
                    if ((i != j)&&(*src >= 0))
                        work_matrix[i * (m.n + 1) + j] = *src;
                    else
                        work_matrix[i * (m.n + 1) + j] = inf;
        }
        //-----------------------------------------------------------------------
        // Начинаем расчёт
        matrix_evaluation(&m, work_matrix, 0, m.n + 1);
        //-----------------------------------------------------------------------
        // Вывод результата
        if (m.global_min < inf)
        {
            p[0] = m.global_min;
            p[1] = m.step;
            unsigned int v_result[n];


            for (unsigned int i = 0; i < n; i++)
            {
                v_result[v_best[i].src - 1] = v_best[i].dest;
            }

            unsigned int j = 0;
            for (unsigned int i = 2; i < m.n + 1; i++)
            {
                p[i] =  v_result[j] - 1;
                j = v_result[j] - 1;
            }
            p[n + 1] = 0;
            return n + 2;
        }
        else
        {
            return 0;
        }
        //-----------------------------------------------------------------------
    }
#ifdef __cplusplus
}
#endif