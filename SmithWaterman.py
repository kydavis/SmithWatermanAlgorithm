
"""
	Use Smith-Waterman algorithm to determine the aligned sequences.

	Parameters:
		seq1: String of first sequence
		seq2: String of second sequence
		sub: Substitution matrix
		gap: gap penalty

	Return:
		Aligned subsequences
"""

class SmithWaterman:

	def __init__(self, sub, gap):
		self.gap = gap
		self.sub = dict(sub)

	def setGap(self, gap):
		self.gap = gap

	def setSubMatrix(self, sub):
		self.sub = sub

	def	_traceback(self, seq1, seq2, scoreMatrix, maxCoord):
		ret1 = []
		ret2 = []
		alignment = []
		i = maxCoord[0]
		j = maxCoord[1]
		score = scoreMatrix[i][j]
		while (score > 0):
			top = scoreMatrix[i - 1][j]
			if (top - score != self.gap):
				top = 0
			left = scoreMatrix[i][j - 1]
			if (left - score != self.gap):
				left = 0
			diag = scoreMatrix[i - 1][j - 1]
			if (score - diag != self.sub[seq1[j - 1]][seq2[i - 1]]):
				diag = 0
			score = max(top, left, diag)
			if (score == diag):
				i = i - 1
				j = j - 1
				ret2.append(seq2[i])
				ret1.append(seq1[j])
				if (seq2[i] == seq1[j]):
					alignment.append("|")
				else:
					alignment.append(":")
			elif (score == top):
				i = i - 1
				ret2.append(seq2[i])
				ret1.append("-")
				alignment.append(" ")
			elif (score == left):
				j = j - 1
				ret1.append(seq1[j])
				ret2.append("-")
				alignment.append(" ")
		ret1 = "".join(ret1[::-1])
		ret2 = "".join(ret2[::-1])
		alignment = "".join(alignment[::-1])
		return ("\n".join([ret1, alignment, ret2]))
			

	def scoreMatrixToStr(self, scoreMatrix):
		ret = []
		for lst in scoreMatrix[0:]:
			ret.append(" ".join(str(x) for x in lst))
		return ("\n".join(ret))

	def	allignSequence(self, seq1, seq2):
		scoreMatrix = [([0] * (len(seq1) + 1)) for row in range(len(seq2) + 1)]
		maxScore = 0
		i = 1
		for n2 in seq2:
			j = 1
			for n1 in seq1:
				gapDown = scoreMatrix[i - 1][j] - self.gap
				gapRight = scoreMatrix[i][j - 1] - self.gap
				match = scoreMatrix[i - 1][j - 1] +  self.sub[n1][n2]
				scoreMatrix[i][j] = max(gapDown, gapRight, match, 0)
				if (scoreMatrix[i][j] > maxScore):
					maxScore = scoreMatrix[i][j]
					maxCoord = [i, j]
				j += 1
			i += 1
		return (self._traceback(seq1, seq2, scoreMatrix, maxCoord))

