#!/usr/bin/env python3

import os
import random

def create_directory(path):
    """Создание директории если она не существует"""
    if not os.path.exists(path):
        os.makedirs(path)

def generate_matrix(filename, rows, cols, min_val=-10.0, max_val=10.0, zero_prob=0.0):
    """Генерация матрицы и запись в файл"""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"{rows} {cols}\n")
        
        for i in range(rows):
            row = []
            for j in range(cols):
                if random.random() < zero_prob:
                    val = 0.0
                else:
                    val = random.uniform(min_val, max_val)
                row.append(f"{val:.6f}")
            f.write(" ".join(row) + "\n")

def main():
    """Основная функция"""
    data_dir = os.path.dirname(os.path.abspath(__file__))
    
    test_cases = [
        # id, rows_a, cols_a, rows_b, cols_b, min_val, max_val, zero_prob
        (1, 2, 3, 3, 2, -5.0, 5.0, 0.0),    # 2x3 * 3x2
        (2, 3, 3, 3, 3, -10.0, 10.0, 0.0),  # 3x3 * 3x3
        (3, 4, 4, 4, 4, -2.0, 2.0, 0.0),    # 4x4 * 4x4
        (4, 2, 2, 2, 2, 0.0, 0.0, 1.0),     # нулевые матрицы
        (5, 2, 2, 2, 2, -5.0, -1.0, 0.0),   # отрицательные значения
        (10, 16, 16, 16, 16, -1.0, 1.0, 0.0),  # для производительности
        (16, 32, 32, 32, 32, -1.0, 1.0, 0.0),  # для производительности
        (20, 64, 64, 64, 64, -1.0, 1.0, 0.0),  # для производительности
        (32, 128, 128, 128, 128, -0.5, 0.5, 0.0),  # для производительности
    ]
    
    for test_id, rows_a, cols_a, rows_b, cols_b, min_val, max_val, zero_prob in test_cases:
        filename_a = os.path.join(data_dir, f"A_{test_id}.txt")
        filename_b = os.path.join(data_dir, f"B_{test_id}.txt")
        
        print(f"Генерация теста {test_id}: {rows_a}x{cols_a} * {rows_b}x{cols_b}")
        
        generate_matrix(filename_a, rows_a, cols_a, min_val, max_val, zero_prob)
        generate_matrix(filename_b, rows_b, cols_b, min_val, max_val, zero_prob)
        
        print(f"  Созданы файлы: {os.path.basename(filename_a)}, {os.path.basename(filename_b)}")
    
    print(f"\nВсе файлы созданы в директории: {data_dir}")

if __name__ == "__main__":
    random.seed(42)  # для воспроизводимости
    main()