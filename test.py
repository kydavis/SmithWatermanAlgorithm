from SmithWaterman import SmithWaterman

subDict = {"A" : {"A" : 3, "T" : -3, "C" : -3, "G" : -3},
			"T": {"A" : -3, "T" : 3, "C" : -3, "G" : -3},
			"C": {"A" : -3, "T" : -3, "C" : 3, "G" : -3},
			"G": {"A" : -3, "T" : -3, "C" : -3, "G" : 3}}

SW = SmithWaterman(subDict, 2)
aligned = SW.allignSequence("TGTTACGG", "GGTTGACTA")
print("Seq1 Raw:{}\nSeq2 Raw:{}".format("TGTTACGG", "GGTTGACTA"))
print(aligned)
