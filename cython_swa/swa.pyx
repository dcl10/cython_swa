import cython
import os

from libc.stdlib cimport malloc, free

cdef int s_ab(str char1, str char2):
    cdef int score = -1
    if char1 == char2:
        score = 1
    return score


def cy_swa(seq1: str, seq2: str):
    """Compute the local alignment of 2 strings using the Smith-Waterman algorithm.

    Args:
        seq1 (str): The first string
        seq2 (str): The second string
    """
    cdef int GAP_PENALTY = 2
    cdef int i = 0, j = 0, current
    cdef Py_ssize_t len_seq1 = len(seq1), len_seq2 = len(seq2)
    cdef int **score_matrix = <int **> malloc((len_seq1 + 1) * (len_seq2 + 1) * sizeof(int))
    for i in range(len_seq1 + 1):
        score_matrix[i] = <int *> malloc((len_seq2 + 1) * sizeof(int))
    
    for i in range(len_seq1 + 1):
        for j in range(len_seq2 + 1):
            score_matrix[i][j] = 0

    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            score_matrix[i][j] = max(
                score_matrix[i - 1][j - 1] + s_ab(seq1[i - 1], seq2[j - 1]),
                score_matrix[i - 1][j] - GAP_PENALTY,
                score_matrix[i][j - 1] - GAP_PENALTY,
                0
            )

    top_str = ""
    mid_str = ""
    btm_str = ""

    cdef Py_ssize_t x = len_seq1, y = len_seq2
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

    free(score_matrix)
    return f"{str.join('', reversed(top_str))}{os.linesep}"\
        f"{str.join('', reversed(mid_str))}{os.linesep}"\
        f"{str.join('', reversed(btm_str))}"
    