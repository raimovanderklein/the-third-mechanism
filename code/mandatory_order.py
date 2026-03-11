#!/usr/bin/env python3
"""
MANDATORY ORDER TEST: Is G→T→A→C the only valid cyclic order?

Tests all 24 possible arrangements of {G, T, A, C} against four
chemical dependency constraints. Only one cyclic order satisfies all four.

From: "The Third Mechanism" (van der Klein, 2026)

Result: Exactly one cycle works: G → T → A → C → G
"""

from itertools import permutations

constraints = [
    ('G', 'T', 'Flexibility required for hydrophobic collapse'),
    ('T', 'A', 'Stable core required for surface function'),
    ('A', 'C', 'Active function required for regulation'),
    ('C', 'G', 'Degradation returns components to free pool'),
]

print("Chemical dependency constraints:")
for before, after, reason in constraints:
    print(f"  {before} → {after}: {reason}")
print()

valid = []
for order in permutations(['G', 'T', 'A', 'C']):
    ok = True
    for before, after, reason in constraints:
        i_b = order.index(before)
        i_a = order.index(after)
        if (i_a - i_b) % 4 != 1:
            ok = False
            break
    if ok:
        valid.append(order)

print(f"Tested {len(list(permutations(['G','T','A','C'])))} possible arrangements.")
print(f"Valid cyclic orders: {len(valid)}")
print()
for v in valid:
    print(f"  {' → '.join(v)} → {v[0]}")
print()

if len(valid) == 4:
    print("All 4 are rotations of the same cycle: G → T → A → C → G")
    print("There is exactly ONE mandatory cyclic order.")
