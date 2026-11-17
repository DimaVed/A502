# TODO: дальнейший порт A502 / PAN AIR

Дата среза: 2026-11-17.

## 1. Что уже готово

- Парсер `$points`, генераторы `$circular` (nac6: 5 сетей, 185 точек,
  128 панелей) и `$trailing` (swb: 9 сетей, 78 панелей), контрольные
  точки.
- Геометрия панелей: P1..P9, 8 подпанелей, collapsed edges.
- `refloc`-фреймы, сжимаемость, `Compressibility`.
- Субпанельные сплайны:
  - source `SPSPLS` (линейный);
  - doublet `SPSPLD` (квадратичный);
  - `alam`, `edgls`.
- Внешние сплайны:
  - `B_S` для source;
  - `B_D` для doublet, включая edge/corner-панели;
  - разреженная сетевая `B_D`.
- PIC:
  - классификация `influ`;
  - subsonic near-field `subsbi` по подпанелям;
  - supersonic subinclined `supsbi` и dispatch `nftpic`;
  - сборка влияния панели и сети;
  - промежуточное `qffcal`-приближение.
- Глобальные AIC-строки, в том числе по нескольким сетям.
- Boundary conditions:
  - ядро `nct=1`, `nct=2`, `nct=4`;
  - локальные source/doublet-члены в произвольной точке панели;
  - маппинг `nlopt` 2..9, 11..20 в коэффициенты BC;
  - правые части `nropt` 0..3, 5..8 (кроме `cbet`/input).
- Matching:
  - полные совпадающие рёбра и общие непрерывные под-цепочки точек;
  - равенство соответствующих doublet-параметров;
  - equivalence classes (union-find), без зависимых строк в общих углах;
  - исключение source-only сетей (`ntd=0`) из doublet-abutments.
- Типы сингулярностей по `kt`, активные неизвестные.
- Таблица default singularity/BC-опций по `kt` (`inbc`).
- Глобальная система:
  - сборка normal-flux BC + matching;
  - выбор активных source/doublet-колонок;
  - row/column equilibration, rank/conditioning-диагностика;
  - аудит числа уравнений/неизвестных по сетям и типам точек;
  - прямой LU-solve с iterative refinement;
  - minimum-norm gauge (COD) для вырожденных/прямоугольных систем;
  - восстановление полного вектора параметров;
  - первый end-to-end `AeroModel::solveFlow()` для одиночной
    source-only/doublet-only сети.
- Постпроцессинг:
  - поверхностные скорости;
  - `cpcal`-формулы: linear/slender/second-order/isentropic;
  - парсинг `$material`/surface assignments и `matvfs`-пропорции
    (сверено с эталонным `nac6.out`);
  - sanity-check на нулевых сингулярностях.

## 2. Известные ограничения

- Для composite-сетей глобальная система недоопределена: не хватает
  части A502-условий и/или исключения зависимых параметров.  Min-norm
  gauge даёт численное решение, но не заменяет эти условия.
- Замкнутая doublet-only сеть имеет постоянный doublet-null-mode
  (velocity-инвариантный); gauge реализован.
- Если навесить normal-flux BC на все 25 контрольных точек замкнутой
  doublet-only сети, система становится несовместной: нужны A502
  edge/corner-условия вместо части этих строк.
- Квадратная doublet-only система плохо обусловлена; физическая невязка
  около `1e-4`.
- Полные `bconcl`/`bcopt` не портированы:
  - есть таблица default-опций по `kt`, но она ещё не применяется при
    сборке строк;
  - нет разбора пользовательских BC-карт и deck-специфичных override
    (`kt=2/18/19`);
  - нет `nlopt` 1, 10, 21..25 (пользовательские коэффициенты, `ccof`)
    и нет `nropt=4` (`cbet`), 9 (input bet);
  - нет `nct=3` и части edge/corner-условий;
  - нет wake-условий (`kt=18/6/20`) в полном объёме.
- Matching поддерживает полные рёбра и непрерывные общие под-цепочки:
  - нет частичных перекрытий с другой дискретизацией/интерполяцией;
  - нет точек пересечения abatment-сегментов;
  - equivalence classes есть для найденных рёбер/углов, но нет полного
    A502-elimination (known/zeroed/dependent параметры);
  - нет полного учёта ориентации/знаков.
- PIC не завершён:
  - нет far-field `ffpic`;
  - нет полного `qnfcal`;
  - `supsbi` (supersonic subinclined) портирован и проверен численным
    интегрированием; нет `supspi`/`aicsup` (superinclined);
  - near-field используется для всех расстояний (медленно и хуже
    обусловлено).
- Материалы парсятся и дают `Cp`-ratios, но ещё не подключены к AIC/BC
  и к выбору формулы `Cp` по сети/поверхности при постпроцессинге.
- `AeroModel::solveFlow()` — только одиночная source-only/doublet-only
  сеть, упрощённый normal-flux BC (`cu=cl=1`, `nropt=2`); composite и
  несколько сетей пока не решаются, `Panel.cp` при отказе остаётся
  заглушкой.
- Из генераторов геометрии есть `$circular` и `$trailing`; нет
  `$quadrilateral`, `$gothic`, `$elliptical`.  В `$circular` не поддержан
  вариант `nopt>0` с z-картами без заголовка `=z`.
- Нет сил/моментов, streamlines, off-body flow field.
- Нет валидации на эталонных результатах PAN AIR.

## 3. Приоритетный план

### Этап 1. Стабилизировать первый solver

Цель: получить устойчивое решение для простых doublet-only и source-only
сетей, затем понять, каких уравнений не хватает composite-сетям.

1. Провести аудит «уравнения vs неизвестные» по типам сетей:
   - для `kt=12/21/22` — doublet-only;
   - для `kt=6/7/17` — source-only;
   - для `kt=1/...` — composite.
2. Применить active-маски непосредственно при сборке AIC/BC-строк,
   а не только при выборе колонок.
3. Добавить масштабирование строк/столбцов и диагностику:
   - число обусловленности;
   - rank;
   - residual;
   - нулевые/зависимые строки.
4. Ввести gauge/regularization для плохо обусловленных квадратных систем.
5. Разобраться, какие дополнительные условия A502 даёт для composite:
   - edge/corner BC;
   - doublet continuity;
   - wake conditions;
   - исключение зависимых параметров.

Критерий: для изолированной `kt=12` системы residual < `1e-8`,
для простой composite-сети матрица становится квадратной или
корректно недоопределённой с явным gauge.

### Этап 2. Полные boundary conditions

1. Портировать обработку опций `bconcl`/`bcopt`:
   - значения по умолчанию для каждого `kt`;
   - пользовательские BC-карты;
   - `nropt=0..9`;
   - `nct=3`;
   - edge/corner-условия.
2. Реализовать RHS от набегающего потока для всех используемых опций.
3. Добавить wake-условия для `kt=18/6/20`.
4. Покрыть тестами каждую поддерживаемую комбинацию `kt`/`nct`/`nropt`.

### Этап 3. Matching и abutments

1. Частичные совпадения рёбер.
2. Несколько сетей в одной точке/угле.
3. Ориентация и знаки doublet-параметров.
4. Elimination зависимых параметров/equivalence classes.
5. Стыковка wake-сетей.
6. Тесты на реальных sample-деках (`btac`, `swb`, `vss`, `ellip`).

### Этап 4. Полный PIC

1. Far-field `ffpic` + необходимые `cf`/`hm`/`hbm`.
2. `qnfcal` для intermediate/quasi-near.
3. Supersonic:
   - `supsbi` — готово;
   - `supspi`/`aicsup` — осталось;
   - `nftpic` dispatch по `iin` — готово.
4. Диспетчеризация по `influenceLevel` вместо near-field везде.
5. Сравнение far/intermediate/near на одних и тех же точках.

### Этап 5. Постпроцессинг и интеграция

1. Second-order и isentropic `Cp` с материалами.
2. Силы и моменты.
3. Streamlines и off-body flow field.
4. Подключить глобальный solve к `AeroModel::solve()`.
5. Записывать реальный `Panel.cp` и результаты в VTK/ResultVisualizer.

### Этап 6. Геометрия и парсер

1. Генераторы `$quadrilateral`, `$gothic`, `$elliptical`
   (`$circular` и `$trailing` уже есть).
2. Полная поддержка `nac6`/соответствующих секций.
3. Недостающие опции `$points`/сетей.

## 4. Ближайшие конкретные задачи

- [x] Применить `hasSource`/`hasDoublet` в
  `evaluateNetworkNearFieldInfluence()`,
  `buildPanelPointLocalTerms()`, `buildGlobalNearFieldAic()`.
- [x] Добавить тесты для `kt=7`, `kt=12` и смешанной системы.
- [x] Добавить row/column scaling и rank/residual-диагностику в
  `GlobalBoundarySystem`.
- [x] Определить и реализовать gauge/regularization для doublet-only.
  Сделано как minimum-norm COD-gauge для вырожденных/прямоугольных систем;
  физический тест на замкнутой коробке выявил несовместность наивного
  набора BC (см. «Известные ограничения»).
- [x] Провести аудит числа уравнений/неизвестных для composite-сетей.
  Вывод: у composite-сети число all-control-point строк равно числу
  doublet-параметров, поэтому дефицит в точности равен числу source-
  параметров; нужен A502-механизм их исключения/доп. условий.
- [ ] Добавить `nct=3` и недостающие edge/corner BC.
- [x] Расширить matching на частичные рёбра.
  Сделаны непрерывные общие под-цепочки точек; произвольные перекрытия
  с интерполяцией остаются в ограничениях.
- [ ] Портировать `ffpic`.
- [ ] Портировать `qnfcal`.
- [x] Портировать `supsbi` и диспетчеризацию `nftpic` по `iin`
  (supersonic subinclined); сверено с численным интегралом ядра.
- [ ] Портировать `supspi`/`aicsup` (supersonic superinclined).
- [x] Подключить solver к `AeroModel`.
  Сделан ограниченный `solveFlow()` (single-network source/doublet-only);
  полный composite-путь ждёт A502 edge/corner-условий и elimination.
- [x] Реализовать second-order/isentropic `Cp`.

## 5. Проверка и команды

Сборка и тесты:

```bash
cmake -S . -B build-test
cmake --build build-test --config Release --target Tests --parallel 4
./build-test/src/Tests/Release/Tests.exe --reporter compact
```

Правила работы:

- один логический шаг — один атомарный коммит;
- новые модули сопровождать тестами;
- не коммитить `build/`, `build-test/`, `logs/`, временные файлы;
- не подменять A502-формулы «правдоподобными» упрощениями без
  явной пометки в коде и тестах.

## 6. Открытые вопросы

- Точный механизм исключения зависимых doublet-параметров в A502.
- Полный набор уравнений для composite-сетей: аудит показал, что при
  строках во всех контрольных точках не хватает ровно по одному условию
  на каждый source-параметр (source и doublet вместе не определяются
  только normal-flux/matching уравнениями).
- Как A502 обрабатывает edge/corner BC и matching в одной точке.
- Параметры материалов и выбор `Cp`-правила по умолчанию.
- Нужна ли поддержка частотного (`dsnic`) режима.
