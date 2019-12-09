// This script contains the construction of 101x101 and 105x105 matrices
//  described in 'Real Lagrangians in Calabi-Yau threefolds'
// (https://arxiv.org/abs/1908.06685) by H. Arguz and T. Prince.
// This source code can be tested at http://magma.maths.usyd.edu.au/calc/
//
// The matrix 'Square' is a 105x105 matrix which stores a matrix representative
// of the map D -> D^2 from H^2(X,Z2) -> H^4(X,Z2), where X is the mirror quintic.
// The vector space H^2(X,Z2) is 101 dimensional, and the triple intersection
// numbers among elements of a spanning set of size 105 in this space 
// are described by Gross in Proposition 4.2 of 'Topological Mirror Symmetry' 
// (https://arxiv.org/abs/math/9909015).

// Gross also specifies triple intersection numbers of a set
// of 101 linearly independent lattice vectors in H^2(X,Z). These do not span H^2(X,Z),
// but descend to a basis of H^2(X,Z2). We store these intersection numbers in Square_101,
// and verify that Square and Square_101 have the same rank.

Z2 :=  FiniteField(2);
Square := Matrix(Z2,105,105,[]);
Square_101 := Matrix(Z2,101,101,[]);

// The 105 cohomology classes L_i, E^l_{i,j}, and E^l_{i,j,k}
// described in 'Topological Mirror Symmetry' correspond to
// the integral points in the boundary of a four-dimensional simplex P.

// We index the vertices of P by elements of {1,2,3,4,5}.
// Edges of P are indexed by pairs (i,j) where i and j are in {1,2,3,4,5} and i < j.
// Triangular faces of P are indexed by a triple (i,j,k) where i, j, and k are in {1,2,3,4,5}, and i < j < k.

// We store the 10 pairs (i,j) indexing the edges of P.
pairs := [[a,b] :  b in [1..5], a in [1..5] | a lt b ];

// We store the 10 triples (i,j,k) indexing the 2-d faces of P.
triples := [[a,b,c] :  c in [1..5], b in [1..5], a in [1..5] | a lt b and b lt c ];

// The four integral points in the relative interior of each edge are indexed by l in {1,2,3,4}.
// These points are labelled in ascending order along the edge from the vertex i to j.

// The six integral points in the relative interior of each triangular face are indexed by l in {1,2,3,4,5,6}.
// These points are indexed as in Figure 4.6 of 'Topological Mirror Symmetry'.
// We note that our indexing of points along the edge between i and k differs from that shown in Figure 4.6:
// Gross orders the points along the edge from vertex k to i in ascending order.
// Under our convention these points are put into ascending order from vertex i to k (since i < k).

// Let L_i denote the divisors corresponding to vertices, indexed by i in {1,2,3,4,5}.
// Let E^l_{i,j}  denote the exceptional divisors corresponding to integral points in edges.
// Let E^l_{i,j,k} denote the exceptional divisors corresponing to integral points in 2-d faces.

// The matrix Square is symmetric, and consists of 9 blocks arranged as shown.
// A  B^t  C^t
// B  D  E^t
// C  E  F
// There are six distinct blocks to construct:
// The matrix A contains the triple intersection numbers (E^{l'}_{i',j'})^2.E^l_{i,j}.
// The matrix B contains the triple intersection numbers (E^{l'}_{i',j'})^2.E^l_{i,j,k}.
// The matrix C contains the triple intersection numbers (E^{l'}_{i',j'})^2.L_a.
// The matrix D contains the triple intersection numbers (E^{l'}_{i',j',k'})^2.E^l_{i,j,k}.
// The matrix E contains the triple intersection numbers (E^{l'}_{i',j',k'})^2.L_a.
// The matrix F contains the triple intersection numbers L_j^2.L_a

// ---------------------------------------------------------------------------------------
// Non-zero entries in the matrix A are the following.
// (1) Entries (E^l_{i,j})^2.E^l_{i,j}, as (E^l_{i,j})^3 = 1 mod 2.
// (2) Entries (E^2_{i,j})^2.E^3_{i,j} and  (E^3_{i,j})^2.E^2_{i,j}.
// (3) Entries (E^l_{i',j'})^2.E^{l'}_{i,j} such that the vertices i, j, i', and j' lie in a common triangle,
// and the integral points corresponding to these divisors are connected by a non-exterior edge in Figure 4.6
// of 'Topological Mirror Symmetry'. Note the difference in our convention for indexing points along edges.
// Under our convention, given a triple i,j,k, the following points are connected by an edge:
// E^1_{i,j} and E^1_{i,k}
// E^4_{i,j} and E^1_{j,k}
// E^4_{j,k} and E^4_{i,k}.

// Entries (1) and (2) correspond to the 4x4 block matrices:
Block := Matrix(Z2,4,4,[[1,0,0,0],[0,1,1,0],[0,1,1,0],[0,0,0,1]]);
for i in [0..9] do
	InsertBlock(~Square,Block,4*i+1,4*i+1);
end for;

// Fixing a triple t = (i,j,k) we find the index of the edges corresponding to pairs
// (t[2],t[3]), (t[1],t[3]), and (t[1],t[2]) in the set pairs. We then assign the non-zero
// entries as described in (3).
for t in triples do
	indices := [Index(pairs,[t[2],t[3]]), Index(pairs,[t[1],t[3]]), Index(pairs,[t[1],t[2]])];

	Square[(indices[2]-1)*4+1,(indices[3]-1)*4+1] := 1;
	Square[(indices[3]-1)*4+1,(indices[2]-1)*4+1] := 1;
	
	Square[(indices[1]-1)*4+1,(indices[3]-1)*4+4] := 1;
	Square[(indices[3]-1)*4+4,(indices[1]-1)*4+1] := 1;
	
	Square[(indices[1]-1)*4+4,(indices[2]-1)*4+4] := 1;
	Square[(indices[2]-1)*4+4,(indices[1]-1)*4+4] := 1;
end for;

print "Check the submatrix A: intersection numbers (E^{l'}_{i',j'})^2.E^l_{i,j}.\n";
print "We check that this submatrix is symmetric:", IsSymmetric(Square);
print "We expect rows E^l_{i,j} to contain 4 non-zero entries if l is in {1,4}. Otherwise we expect rows E^l_{i,j} to contain 2 non-zero entries.";
number_of_entries := [[ #[1: k in [1..40] | Square[i*4+j,k] eq 1] : j in [1..4]] : i in [0..9]];
print number_of_entries eq [[4,2,2,4] : i in [0..9]];
print "----------\n";

// ---------------------------------------------------------------------------------------
// The submatrix B is a 60x40 matrix divided into 100 6x4 blocks. Each of these subblocks
// is one of four possible matrices: the zero matrix, or one of the three blocks given below.
// These are read from Figure 4.6 of 'Topological Mirror Symmetry'. For example, given a triple
// (i,j,k) and pair (j,k) the points indexed by 1,2,3,4 along the edge (j,k) are connected by
// a non-exterior edge to {4}, {4,5}, {5,6} and {6} respectively. Hence the corresponding block matrix stores
// non-zero entries in positions (4,1), (4,2), (5,2), (5,3), (6,3), and (6,4).

Block := Matrix(Z2,60,40,[]);
MiniBlocks := [Matrix(Z2,6,4,[[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,1,0,0],[0,1,1,0],[0,0,1,1]]),
Matrix(Z2,6,4,[[1,1,0,0],[0,0,0,0],[0,1,1,0],[0,0,0,0],[0,0,0,0],[0,0,1,1]]),
Matrix(Z2,6,4,[[1,1,0,0],[0,1,1,0],[0,0,0,0],[0,0,1,1],[0,0,0,0],[0,0,0,0]])
];

// r_pos stores the index in [1..60] of the top row of the current subblock.
r_pos := 1;

for t in triples do

	// c_pos stores the index in [1..40] of the left most column of the current subblock.
	c_pos := 1;

		for p in pairs do

		// If the pair p = [a,b] is contained in the triple t = [c,d,e], there is a non-zero
		// block with top left position (r_pos,c_pos). We determine which block to add
		// by computing which index of the triple [c,d,e] which does not appear in the pair [a,b].
		if SequenceToSet(p) subset SequenceToSet(t) then
			mini_block :=  [k : k in [1,2,3] | t[k] notin p][1];
			InsertBlock(~Block,MiniBlocks[mini_block],r_pos,c_pos);
		end if;

		c_pos := c_pos + 4;
	end for;
	r_pos := r_pos + 6;
end for;


print "Check the submatrix B: intersection numbers (E^{l'}_{i',j'})^2.E^l_{i,j,k}.\n";
print "We expect row E^l_{i,j,k} to contain 4 non-zero entries if l is in {1,4,6}. Otherwise we expect row E^l_{i,j,k} to contain 2 non-zero entries.";
number_of_entries := [[ #[1: k in [1..40] | Block[i*6+j,k] eq 1] : j in [1..6]] : i in [0..9]];
print number_of_entries eq [[4,2,2,4,2,4] : i in [0..9]];
print "We expect column E^l_{i,j} to contain 6 non-zero entries if l is in {2,3}. Otherwise we expect column E^l_{i,j} to contain 3 non-zero entries.";
number_of_entries := [[ #[1: k in [1..60] | Block[k,i*4+j] eq 1] : j in [1..4]] : i in [0..9]];
print number_of_entries eq [[3,6,6,3] : i in [0..9]];
print "----------\n";

InsertBlock(~Square,Block,41,1);
InsertBlock(~Square,Transpose(Block),1,41);

// ---------------------------------------------------------------------------------------
// Insert the matrix blocks C and C^t.
// Nonzero entries are given by the intersection numbers L_i^2.E^1_{i,j} = 1.
// Nonzero entries are given by the intersection numbers L_j^2.E^4_{i,j} = 1.
// Note that, the convention in 'Topological Mirror Symmetry' is to identify 
// L_j^2.E^4_{i,j} = L_j^2.E^1_{j,i}, and records that L_j^2.E^1_{j,i} = 1 mod 2.

Block := Matrix(Z2,5,40,[]);

for i in [1..5] do
	
	// c_pos stores the index of the left most column of each block.
	c_pos := 1;
	
	// Each of the 10 blocks which comprise 'Block'
	// is indexed by a pair [a,b], where a < b.
	for p in pairs do

			// If i is equal to a, then record that L_i^2.E^1_{i,j} = 1.
			if i eq p[1] then
				Block[i,c_pos] := 1;
			end if;

			// If i is equal to b, then record that L_j^2.E^4_{i,j} = 1.
			if i eq p[2] then
				Block[i,c_pos+3] := 1;
			end if;
			c_pos := c_pos + 4;
	end for;
end for;

print "Check the submatrix C: intersection numbers (E^{l'}_{i',j'})^2.L_a.\n";
print "We expect each row L_i to contain precisely 4 non-zero entries";
print [ #[1 : j in [1..40] | Block[i,j] eq 1] : i in [1..5]] eq [4,4,4,4,4],"\n";
print "We expect column E^l_{i,j} to have a non-zero entry if and only if l is 1 or 4.";
number_of_entries := [ [ #[1 : k in [1..5] | Block[k,i*4+j] eq 1] : j in [1..4] ] : i in [0..9]];
print number_of_entries eq [[1,0,0,1] : i in [0..9]];
print "----------\n";

InsertBlock(~Square,Block,101,1);
InsertBlock(~Square,Transpose(Block),1,101);

// ---------------------------------------------------------------------------------------
// Insert the matrix block D.
// This matrix contains a block for each triple (i,j,k).
// We read this diectly from Figure 4.6 of 'Topological Mirror Symmetry'.
Block := Matrix(Z2,6,6,[[0,1,1,0,0,0],[1,0,1,1,1,0],[1,1,0,0,1,1],[0,1,0,0,1,0],[0,1,1,1,0,1],[0,0,1,0,1,0]]);
for i in [0..9] do
	InsertBlock(~Square,Block,6*i+41,6*i+41);
end for;

print "Check the submatrix D: intersection numbers (E^{l'}_{i',j',k'})^2.E^l_{i,j,k}.\n";
print "We expect column E^l_{i,j,k} to contain 2 non-zero entries if l is 1, 4, or 6. Otherwise we expect column E^l_{i,j,k} to contain 4 non-zero entries.";
number_of_entries :=  [[ #[1 : k in [1..6] | Square[i*6+j+40,i*6+k+40] eq 1 ] : j in [1..6]] : i in [0..9]];
print number_of_entries eq [ [2,4,4,2,4,2] : i in [0..9]];
print "----------\n";

// ---------------------------------------------------------------------------------------

// The matrix blocks E and E^t contain only zero entries.
print "Check the submatrix E: intersection numbers (E^{l'}_{i',j',k'})^2.L_a.\n";
print "We expect blocks E and its transpose to contain only zero entries.";
print ZeroMatrix(Z2,5,60) eq Submatrix(Square,101,41,5,60) and ZeroMatrix(Z2,60,5) eq Submatrix(Square,41,101,60,5);
print "----------\n";

// ---------------------------------------------------------------------------------------
// The matrix F is a square identity matrix since L_i^3 is non-zero mod 2,
// and L^2_iL_j = 0 mod 2 for any distinct i and j.

InsertBlock(~Square,IdentityMatrix(Z2,5),101,101);

print "Check the submatrix F: intersection numbers (L_j)^2.L_a.\n";
print "We expect this to be a 5x5 identity block.";
print IdentityMatrix(Z2,5) eq Submatrix(Square,101,101,5,5);
print "----------\n";

print "Check Square is a symmetric matrix:", IsSymmetric(Square), "\n";
// We compute the rank of this matrix.
print "Rank of Square:", Rank(Square);

// ---------------------------------------------------------------------------------------

// We now construct the matrix Square_101.
// This matrix consists of triple intersection numbers between the classes
// E^l_{i,j}, E^l_{i,j,k} and H, as described in Proposition 4.2 of 'Topological Mirror Symmetry'.
// Since the first 100 rows/columns of the matrix 'Square' store the triple intersection numbers 
// between E^l_{i,j} and E^l_{i,j,k} these entries are unchanged in Square_101.

InsertBlock(~Square_101, Submatrix(Square,1,1,100,100),1,1);

// By Proposition 4.2 of 'Topological Mirror Symmetry' H^3 = 1 mod 2, and H^2.D = D^2.H = 0 mod 2
// where D is any class E^l_{i,j} or E^l_{i,j,k}. Hence the only non-zero element in the 101st row/column
// is the entry at position (101,101), corresponding to H^2.H.

Square_101[101,101] := 1;
print "\nCheck Square_101 is a symmetric matrix:", IsSymmetric(Square_101);
print "Rank of Square_101:", Rank(Square_101);
print "Matrices Square and Square_101 have equal rank:", Rank(Square) eq Rank(Square_101);
