import sys

def load_cnf(clauses, file_name):
	f = open(file_name, 'r')
	if not f:
		print("Failed to open file {}".format(file_name))
		exit()
	while True:
		line = f.readline().strip()
		if not line: break
		if line[0] == 'p' or line[0] == 'c': continue
		clause = [int(l) for l in line.split()[:-1]]
		clause.sort()
		clauses.append(clause)
	f.close()

# Find longest intersection between two clauses
def intersect(c1, c2):
	i1 = 0
	i2 = 0
	intersection = []
	while i1 < len(c1) and i2 < len(c2):
		if c1[i1] < c2[i2]:
			i1 += 1
		elif c1[i1] > c2[i2]:
			i2 += 1
		else:
			intersection.append(c1[i1])
			i1 += 1
			i2 += 1
	return intersection

# Count the number of times each subexpression occurs
def count_subexprs(subexprs, clauses, min_subexpr_len = 1):
	for i1, c1 in enumerate(clauses):
		found = set()
		for i2, c2 in enumerate(clauses):
			if i1 == i2: continue
			intersection = frozenset(intersect(c1, c2))
			if (
				len(intersection) >= min_subexpr_len and
				intersection not in found
			):
				found.add(intersection)
				if intersection in subexprs:
					subexprs[intersection] += 1
				else: subexprs[intersection] = 1
	return subexprs

def main(base_cnf, learnt_cnf, min_subexpr_len):
	# Load clauses
	clauses = []
	load_cnf(clauses, base_cnf)
	load_cnf(clauses, learnt_cnf)

	# Count subexpressions
	subexprs = {}
	count_subexprs(subexprs, clauses, min_subexpr_len)

	# Sort counts by count
	subexprs = {k: v for k, v in sorted(subexprs.items(), key=lambda item: -item[1])}

	# Output counts
	print("Found {} distinct subexpressions with length >= {}".format(len(subexprs), min_subexpr_len))
	num_to_output = 10
	for i, item in enumerate(subexprs.items()):
		if i == num_to_output: break
		keystr = ' '.join([str(j) for j in item[0]])
		print("{}: {}".format(keystr, item[1]))

if __name__ == '__main__':
	base_cnf   = sys.argv[1]
	learnt_cnf = base_cnf[0:-len('.cnf')] + '_learnt.cnf'
	min_subexpr_len = int(sys.argv[2])

	main(base_cnf, learnt_cnf, min_subexpr_len)