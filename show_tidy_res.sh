#!/bin/bash
# остановка при ошибке любой команды
set -e

# Проверка наличия clang-tidy-21
if ! command -v clang-tidy-21 &> /dev/null; then
    echo "Ошибка: clang-tidy-21 не найден. Установите его:"
    echo "  sudo apt-get install -y clang-tidy-21"
    exit 1
fi

# Проверка наличия run-clang-tidy-21
if ! command -v run-clang-tidy-21 &> /dev/null; then
    echo "Ошибка: run-clang-tidy-21 не найден. Установите его:"
    echo "  sudo apt-get install -y clang-tidy-21"
    exit 1
fi

echo "Используется clang-tidy-21 версии: $(clang-tidy-21 --version | head -n1)"

# Сборка проекта
echo "Сборка проекта..."
cmake --build build --parallel 1

# Создаем директорию для логов, если её нет
mkdir -p tidy_logs

# Удаляем старый tidy.log, если есть
rm -f tidy.log

# Запуск clang-tidy-21
echo "Запуск clang-tidy-21..."
run-clang-tidy-21 -p build > tidy.log 2>&1 || true

# Удаляем старый mytidyprob.log, если есть
rm -f mytidyprob.log

# Фильтруем строки с "liulin"
if [ -s tidy.log ]; then
    grep "liulin" tidy.log > mytidyprob.log || touch mytidyprob.log
else
    touch mytidyprob.log
fi

# Разделяем логи по категориям
echo "Разделение логов по категориям..."

rm -f ./tidy_logs/test_func.log
grep "functional" mytidyprob.log > ./tidy_logs/test_func.log || touch ./tidy_logs/test_func.log

rm -f ./tidy_logs/test_perf.log
grep "performance" mytidyprob.log > ./tidy_logs/test_perf.log || touch ./tidy_logs/test_perf.log

rm -f ./tidy_logs/test_mpi.log
grep "mpi" mytidyprob.log > ./tidy_logs/test_mpi.log || touch ./tidy_logs/test_mpi.log

rm -f ./tidy_logs/test_seq.log
grep "seq" mytidyprob.log > ./tidy_logs/test_seq.log || touch ./tidy_logs/test_seq.log

rm -f ./tidy_logs/test_common.log
grep "common" mytidyprob.log > ./tidy_logs/test_common.log || touch ./tidy_logs/test_common.log

# Подсчет проблем
TOTAL_ISSUES=$(grep -c "error:" mytidyprob.log 2>/dev/null || echo "0")
FUNC_ISSUES=$(grep -c "error:" ./tidy_logs/test_func.log 2>/dev/null || echo "0")
PERF_ISSUES=$(grep -c "error:" ./tidy_logs/test_perf.log 2>/dev/null || echo "0")
MPI_ISSUES=$(grep -c "error:" ./tidy_logs/test_mpi.log 2>/dev/null || echo "0")
SEQ_ISSUES=$(grep -c "error:" ./tidy_logs/test_seq.log 2>/dev/null || echo "0")
COMMON_ISSUES=$(grep -c "error:" ./tidy_logs/test_common.log 2>/dev/null || echo "0")

# Очистка временных файлов
rm -f mytidyprob.log tidy.log

echo ""
echo "=========================================="
echo "Clang-tidy-21 завершен!"
echo "=========================================="
echo "Всего проблем с 'liulin': $TOTAL_ISSUES"
echo "  - functional: $FUNC_ISSUES"
echo "  - performance: $PERF_ISSUES"
echo "  - mpi: $MPI_ISSUES"
echo "  - seq: $SEQ_ISSUES"
echo "  - common: $COMMON_ISSUES"
echo "=========================================="
echo "Логи сохранены в tidy_logs/"

