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

class Counter:
	def __init__(self, clause_id):
		self.count         = 0
		self.min_clause_id = clause_id
		self.max_clause_id = clause_id

	def update(self, clause_id):
		self.count += 1
		self.min_clause_id = min(self.min_clause_id, clause_id)
		self.max_clause_id = max(self.max_clause_id, clause_id)

	def span(self):
		return self.max_clause_id - self.min_clause_id

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
				if intersection not in subexprs: subexprs[intersection] = Counter(i1)
				subexprs[intersection].update(i2)
	return subexprs

def main(base_cnf, learnt_cnf, min_subexpr_len):
	# Load clauses
	clauses = []
	load_cnf(clauses, base_cnf)
	load_cnf(clauses, learnt_cnf)

	# Count subexpressions
	subexprs = {}
	count_subexprs(subexprs, clauses, min_subexpr_len)

	# Calculate totals
	maxLen = max(len(key) for key in subexprs.keys())
	totals = []
	for i in range(maxLen): totals.append([0, 0., 0.]) # num, count, span
	for item in subexprs.items():
		totals[len(item[0]) - 1][0] += 1
		totals[len(item[0]) - 1][1] += item[1].count
		totals[len(item[0]) - 1][2] += item[1].span()
	
	# Output averages
	print("Found {} distinct subexpressions with length >= {}".format(len(subexprs), min_subexpr_len))
	print("length, num, avg_count, avg_span")
	for i, total in enumerate(totals):
		if total[0] == 0: print("{} {} {} {}".format(i + 1, total[0], total[1], total[2]))
		else:             print("{} {} {} {}".format(i + 1, total[0], total[1] / total[0], total[2] / total[0]))

if __name__ == '__main__':
	base_cnf        = sys.argv[1]
	learnt_cnf      = base_cnf[0:-len('.cnf')] + '_learnt.cnf'
	min_subexpr_len = int(sys.argv[2]) #  1 means output everything

	main(base_cnf, learnt_cnf, min_subexpr_len)