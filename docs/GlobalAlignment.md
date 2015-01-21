##Approach

Pretty straight forward.
Assuming M[i][j] represents the optimal score for 
aligning two sequences. Then M[i][j]
is given by the famous Needleman Wunsch Dynamic programming
algorithm:

M[i][j] = max(M_mm, M__i1, M_i2)

where M_mm = M[i-1][j-1] + W(U_i, U_j)
M_i1 = M[i-1][j] + W(U_i, -) => Aligning U(U_i) against a gap
M_i1 = M[i-1][j] + W(U-, V_j) =>  Aligning V(V_j) against a gap

## Initialization

We add “_ _ ” before each sequence
nad 
