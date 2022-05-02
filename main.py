import os
import time

from cython_swa.swa import cy_swa


def py_swa(seq1: str, seq2: str) -> None:
    """Compute the local alignment of 2 strings using the Smith-Waterman algorithm.

    Args:
        seq1 (str): The first string
        seq2 (str): The second string
    """
    GAP_PENALTY = 2
    
    score_matrix = [
        [0 for j in range(len(seq2) + 1)] 
        for i in range(len(seq1) + 1)
    ]

    s_ab = lambda a, b: 1 if a == b else -1
    
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score_matrix[i][j] = max(
                score_matrix[i - 1][j - 1] + s_ab(seq1[i - 1], seq2[j - 1]),
                score_matrix[i - 1][j] - GAP_PENALTY,
                score_matrix[i][j - 1] - GAP_PENALTY,
                0
            )

    top_str = ""
    mid_str = ""
    btm_str = ""

    x, y = len(seq1), len(seq2)
    while x > 0 and y > 0:
        score = max(
            score_matrix[x - 1][y - 1],
            score_matrix[x - 1][y],
            score_matrix[x][y - 1],
        )
        if score == score_matrix[x - 1][y - 1]:
            top_str += seq1[x - 1]
            mid_str += "|"
            btm_str += seq2[y - 1]
            x -= 1
            y -= 1
        elif score == score_matrix[x - 1][y]:
            top_str += seq1[x - 1]
            mid_str += " "
            btm_str += "-"
            x -= 1
        elif score == score_matrix[x][y - 1]:
            top_str += "-"
            mid_str += " "
            btm_str += seq2[y - 1]
            y -= 1

    return f"{str.join('', reversed(top_str))}{os.linesep}"\
        f"{str.join('', reversed(mid_str))}{os.linesep}"\
        f"{str.join('', reversed(btm_str))}"


if __name__ == "__main__":
    start = time.process_time()
    res = py_swa("AAATTTA"*1000, "AAATT"*1000)
    end = time.process_time()
    py_time = end - start
    print(f"Python SWA took {(py_time):.6f}s")

    start = time.process_time()
    res = cy_swa("AAATTTA"*1000, "AAATT"*1000)
    end = time.process_time()
    cy_time = end - start
    print(f"Cython SWA took {(cy_time):.6f}s")

    print(f"Cython speed up of {(py_time / cy_time):.2f}x")
