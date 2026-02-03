"""Simple Needleman-Wunsch sequence alignment implementation.

This module provides global sequence alignment using the classic
Needleman-Wunsch dynamic programming algorithm.
"""

from typing import Tuple, List


# BLOSUM62 scoring matrix (simplified)
BLOSUM62 = {
    ('A', 'A'): 4,  ('A', 'C'): 0,  ('A', 'D'): -2, ('A', 'E'): -1, ('A', 'F'): -2,
    ('A', 'G'): 0,  ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'K'): -1, ('A', 'L'): -1,
    ('A', 'M'): -1, ('A', 'N'): -2, ('A', 'P'): -1, ('A', 'Q'): -1, ('A', 'R'): -1,
    ('A', 'S'): 1,  ('A', 'T'): 0,  ('A', 'V'): 0,  ('A', 'W'): -3, ('A', 'Y'): -2,

    ('C', 'C'): 9,  ('C', 'D'): -3, ('C', 'E'): -4, ('C', 'F'): -2, ('C', 'G'): -3,
    ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'K'): -3, ('C', 'L'): -1, ('C', 'M'): -1,
    ('C', 'N'): -3, ('C', 'P'): -3, ('C', 'Q'): -3, ('C', 'R'): -3, ('C', 'S'): -1,
    ('C', 'T'): -1, ('C', 'V'): -1, ('C', 'W'): -2, ('C', 'Y'): -2,

    ('D', 'D'): 6,  ('D', 'E'): 2,  ('D', 'F'): -3, ('D', 'G'): -1, ('D', 'H'): -1,
    ('D', 'I'): -3, ('D', 'K'): -1, ('D', 'L'): -4, ('D', 'M'): -3, ('D', 'N'): 1,
    ('D', 'P'): -1, ('D', 'Q'): 0,  ('D', 'R'): -2, ('D', 'S'): 0,  ('D', 'T'): -1,
    ('D', 'V'): -3, ('D', 'W'): -4, ('D', 'Y'): -3,

    ('E', 'E'): 5,  ('E', 'F'): -3, ('E', 'G'): -2, ('E', 'H'): 0,  ('E', 'I'): -3,
    ('E', 'K'): 1,  ('E', 'L'): -3, ('E', 'M'): -2, ('E', 'N'): 0,  ('E', 'P'): -1,
    ('E', 'Q'): 2,  ('E', 'R'): 0,  ('E', 'S'): 0,  ('E', 'T'): -1, ('E', 'V'): -2,
    ('E', 'W'): -3, ('E', 'Y'): -2,

    ('F', 'F'): 6,  ('F', 'G'): -3, ('F', 'H'): -1, ('F', 'I'): 0,  ('F', 'K'): -3,
    ('F', 'L'): 0,  ('F', 'M'): 0,  ('F', 'N'): -3, ('F', 'P'): -4, ('F', 'Q'): -3,
    ('F', 'R'): -3, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'V'): -1, ('F', 'W'): 1,
    ('F', 'Y'): 3,

    ('G', 'G'): 6,  ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'K'): -2, ('G', 'L'): -4,
    ('G', 'M'): -3, ('G', 'N'): 0,  ('G', 'P'): -2, ('G', 'Q'): -2, ('G', 'R'): -2,
    ('G', 'S'): 0,  ('G', 'T'): -2, ('G', 'V'): -3, ('G', 'W'): -2, ('G', 'Y'): -3,

    ('H', 'H'): 8,  ('H', 'I'): -3, ('H', 'K'): -1, ('H', 'L'): -3, ('H', 'M'): -2,
    ('H', 'N'): 1,  ('H', 'P'): -2, ('H', 'Q'): 0,  ('H', 'R'): 0,  ('H', 'S'): -1,
    ('H', 'T'): -2, ('H', 'V'): -3, ('H', 'W'): -2, ('H', 'Y'): 2,

    ('I', 'I'): 4,  ('I', 'K'): -3, ('I', 'L'): 2,  ('I', 'M'): 1,  ('I', 'N'): -3,
    ('I', 'P'): -3, ('I', 'Q'): -3, ('I', 'R'): -3, ('I', 'S'): -2, ('I', 'T'): -1,
    ('I', 'V'): 3,  ('I', 'W'): -3, ('I', 'Y'): -1,

    ('K', 'K'): 5,  ('K', 'L'): -2, ('K', 'M'): -1, ('K', 'N'): 0,  ('K', 'P'): -1,
    ('K', 'Q'): 1,  ('K', 'R'): 2,  ('K', 'S'): 0,  ('K', 'T'): -1, ('K', 'V'): -2,
    ('K', 'W'): -3, ('K', 'Y'): -2,

    ('L', 'L'): 4,  ('L', 'M'): 2,  ('L', 'N'): -3, ('L', 'P'): -3, ('L', 'Q'): -2,
    ('L', 'R'): -2, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'V'): 1,  ('L', 'W'): -2,
    ('L', 'Y'): -1,

    ('M', 'M'): 5,  ('M', 'N'): -2, ('M', 'P'): -2, ('M', 'Q'): 0,  ('M', 'R'): -1,
    ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'V'): 1,  ('M', 'W'): -1, ('M', 'Y'): -1,

    ('N', 'N'): 6,  ('N', 'P'): -2, ('N', 'Q'): 0,  ('N', 'R'): 0,  ('N', 'S'): 1,
    ('N', 'T'): 0,  ('N', 'V'): -3, ('N', 'W'): -4, ('N', 'Y'): -2,

    ('P', 'P'): 7,  ('P', 'Q'): -1, ('P', 'R'): -2, ('P', 'S'): -1, ('P', 'T'): -1,
    ('P', 'V'): -2, ('P', 'W'): -4, ('P', 'Y'): -3,

    ('Q', 'Q'): 5,  ('Q', 'R'): 1,  ('Q', 'S'): 0,  ('Q', 'T'): -1, ('Q', 'V'): -2,
    ('Q', 'W'): -2, ('Q', 'Y'): -1,

    ('R', 'R'): 5,  ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'V'): -3, ('R', 'W'): -3,
    ('R', 'Y'): -2,

    ('S', 'S'): 4,  ('S', 'T'): 1,  ('S', 'V'): -2, ('S', 'W'): -3, ('S', 'Y'): -2,

    ('T', 'T'): 5,  ('T', 'V'): 0,  ('T', 'W'): -2, ('T', 'Y'): -2,

    ('V', 'V'): 4,  ('V', 'W'): -3, ('V', 'Y'): -1,

    ('W', 'W'): 11, ('W', 'Y'): 2,

    ('Y', 'Y'): 7,
}

# Make matrix symmetric
for (aa1, aa2), score in list(BLOSUM62.items()):
    if (aa2, aa1) not in BLOSUM62:
        BLOSUM62[(aa2, aa1)] = score


def get_score(aa1: str, aa2: str) -> int:
    """Get BLOSUM62 score for two amino acids."""
    aa1 = aa1.upper()
    aa2 = aa2.upper()
    return BLOSUM62.get((aa1, aa2), -4)  # Default penalty for unknown pairs


def needleman_wunsch(
    seq1: str,
    seq2: str,
    gap_open: int = -10,
    gap_extend: int = -1,
    use_affine_gap: bool = True
) -> Tuple[str, str, float]:
    """
    Perform global sequence alignment using Needleman-Wunsch algorithm.

    Args:
        seq1: First sequence
        seq2: Second sequence
        gap_open: Gap opening penalty (negative value)
        gap_extend: Gap extension penalty (negative value)
        use_affine_gap: Use affine gap penalty (separate open/extend costs)

    Returns:
        Tuple of (aligned_seq1, aligned_seq2, score)
    """
    m, n = len(seq1), len(seq2)

    if use_affine_gap:
        # Three matrices for affine gap penalty
        # M: match/mismatch, Ix: gap in seq1, Iy: gap in seq2
        M = [[float('-inf')] * (n + 1) for _ in range(m + 1)]
        Ix = [[float('-inf')] * (n + 1) for _ in range(m + 1)]
        Iy = [[float('-inf')] * (n + 1) for _ in range(m + 1)]

        # Traceback matrices
        trace_M = [[None] * (n + 1) for _ in range(m + 1)]
        trace_Ix = [[None] * (n + 1) for _ in range(m + 1)]
        trace_Iy = [[None] * (n + 1) for _ in range(m + 1)]

        # Initialize
        M[0][0] = 0
        for i in range(1, m + 1):
            Iy[i][0] = gap_open + (i - 1) * gap_extend
            trace_Iy[i][0] = ('Iy', i - 1, 0)
        for j in range(1, n + 1):
            Ix[0][j] = gap_open + (j - 1) * gap_extend
            trace_Ix[0][j] = ('Ix', 0, j - 1)

        # Fill matrices
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # M matrix (match/mismatch)
                match_score = get_score(seq1[i-1], seq2[j-1])
                from_M = M[i-1][j-1] + match_score
                from_Ix = Ix[i-1][j-1] + match_score
                from_Iy = Iy[i-1][j-1] + match_score

                M[i][j] = max(from_M, from_Ix, from_Iy)
                if M[i][j] == from_M:
                    trace_M[i][j] = ('M', i - 1, j - 1)
                elif M[i][j] == from_Ix:
                    trace_M[i][j] = ('Ix', i - 1, j - 1)
                else:
                    trace_M[i][j] = ('Iy', i - 1, j - 1)

                # Ix matrix (gap in seq1)
                from_M_open = M[i][j-1] + gap_open
                from_Ix_extend = Ix[i][j-1] + gap_extend

                Ix[i][j] = max(from_M_open, from_Ix_extend)
                if Ix[i][j] == from_M_open:
                    trace_Ix[i][j] = ('M', i, j - 1)
                else:
                    trace_Ix[i][j] = ('Ix', i, j - 1)

                # Iy matrix (gap in seq2)
                from_M_open = M[i-1][j] + gap_open
                from_Iy_extend = Iy[i-1][j] + gap_extend

                Iy[i][j] = max(from_M_open, from_Iy_extend)
                if Iy[i][j] == from_M_open:
                    trace_Iy[i][j] = ('M', i - 1, j)
                else:
                    trace_Iy[i][j] = ('Iy', i - 1, j)

        # Find best final state
        final_score = max(M[m][n], Ix[m][n], Iy[m][n])
        if final_score == M[m][n]:
            current_matrix = 'M'
            trace = trace_M
        elif final_score == Ix[m][n]:
            current_matrix = 'Ix'
            trace = trace_Ix
        else:
            current_matrix = 'Iy'
            trace = trace_Iy

        # Traceback
        aln1, aln2 = [], []
        i, j = m, n

        while i > 0 or j > 0:
            if current_matrix == 'M':
                aln1.append(seq1[i-1])
                aln2.append(seq2[j-1])
                if trace[i][j]:
                    current_matrix, i, j = trace[i][j]
                else:
                    break
                trace = trace_M if current_matrix == 'M' else (trace_Ix if current_matrix == 'Ix' else trace_Iy)
            elif current_matrix == 'Ix':
                aln1.append('-')
                aln2.append(seq2[j-1])
                if trace_Ix[i][j]:
                    current_matrix, i, j = trace_Ix[i][j]
                else:
                    break
                trace = trace_M if current_matrix == 'M' else trace_Ix
            else:  # Iy
                aln1.append(seq1[i-1])
                aln2.append('-')
                if trace_Iy[i][j]:
                    current_matrix, i, j = trace_Iy[i][j]
                else:
                    break
                trace = trace_M if current_matrix == 'M' else trace_Iy

        return ''.join(reversed(aln1)), ''.join(reversed(aln2)), final_score

    else:
        # Simple linear gap penalty
        dp = [[0] * (n + 1) for _ in range(m + 1)]
        trace = [[None] * (n + 1) for _ in range(m + 1)]

        # Initialize first row and column
        for i in range(m + 1):
            dp[i][0] = gap_extend * i
        for j in range(n + 1):
            dp[0][j] = gap_extend * j

        # Fill DP table
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = dp[i-1][j-1] + get_score(seq1[i-1], seq2[j-1])
                delete = dp[i-1][j] + gap_extend
                insert = dp[i][j-1] + gap_extend

                dp[i][j] = max(match, delete, insert)

                if dp[i][j] == match:
                    trace[i][j] = 'M'
                elif dp[i][j] == delete:
                    trace[i][j] = 'D'
                else:
                    trace[i][j] = 'I'

        # Traceback
        aln1, aln2 = [], []
        i, j = m, n

        while i > 0 or j > 0:
            if i > 0 and j > 0 and trace[i][j] == 'M':
                aln1.append(seq1[i-1])
                aln2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif i > 0 and trace[i][j] == 'D':
                aln1.append(seq1[i-1])
                aln2.append('-')
                i -= 1
            else:
                aln1.append('-')
                aln2.append(seq2[j-1])
                j -= 1

        return ''.join(reversed(aln1)), ''.join(reversed(aln2)), dp[m][n]


def compute_identity(aln1: str, aln2: str) -> float:
    """Calculate percent identity from aligned sequences."""
    if len(aln1) != len(aln2):
        raise ValueError("Aligned sequences must have same length")

    matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != '-')
    aligned_len = sum(1 for a, b in zip(aln1, aln2) if a != '-' or b != '-')

    return matches / aligned_len if aligned_len > 0 else 0.0


def alignment_to_columns(aln1: str, aln2: str) -> List[Tuple[int, int, int]]:
    """
    Convert aligned sequences to column format like foldseek.

    Returns:
        List of (col, qpos, tpos) tuples where positions are 1-indexed
        or None if gap
    """
    cols = []
    qpos = 0
    tpos = 0

    for col_idx, (aa1, aa2) in enumerate(zip(aln1, aln2), 1):
        q = None if aa1 == '-' else (qpos + 1)
        t = None if aa2 == '-' else (tpos + 1)

        cols.append((col_idx, q, t))

        if aa1 != '-':
            qpos += 1
        if aa2 != '-':
            tpos += 1

    return cols
