from Bio import SeqIO


templates = list(SeqIO.parse("amplicons.fasta", "fasta"))
pre_templates = [str(t.seq) for t in templates]
template_id = [str(t.id) for t in templates]

print(pre_templates[1])

# k = 77

# min_score = int(k)


# kmers = [pre_templates[0][i:i + k] for i in range(0, len(pre_templates[0]) - k + 1)]

# scores = []
file = open("found_templates.fasta","w+")
for index,t in enumerate(pre_templates):
    file.write(">" + (template_id[index])[:-6] + "\n")


    # print((template_id[index])[:-6])

    candidates = {}
    for k in range(3, len(pre_templates[index])):
        min_score = int(k * 0.75)
        kmers = [pre_templates[index][i:i + k] for i in range(0, len(pre_templates[index]) - k + 1)]
        scores = []
        repeats = []

        for kmer_index in range(len(kmers) - k):
            listOfSplits = [pre_templates[index][i:i + k] for i in range(kmer_index, len(pre_templates[index]), k)]


            #Get rid on incomplete splits
            if len(listOfSplits[-1]) != len(listOfSplits[0]):
                listOfSplits = listOfSplits[0:len(listOfSplits)-1]

            #Get rid of first split(it's the same as the k-mer)
            listOfSplits = listOfSplits[1:]


            idx = 0

            kmer_score = 0
            number_of_repeats = 0

            for split in listOfSplits:
                # print(idx)
                idx += 1
                score = sum(el1 == el2 for el1, el2 in zip(kmers[kmer_index], split))
                # print(score)

                if score < min_score:
                    # print("broke")
                    break

                if score >= min_score:
                    kmer_score += score
                    number_of_repeats += 1
                    # print(kmer_score)
                    # print(kmers[kmer_index])

            scores.append(kmer_score)
            repeats.append(number_of_repeats)

        if len(scores) == 0:
            max_score = 0

        if len(scores) > 0:
            max_score = max(scores)

            if max_score != 0:
                idx_max_score = scores.index(max_score)
                candidates[kmers[idx_max_score]] = max_score


        # print(repeats)
        # print(max_score)

        # print(kmers[idx_max_score])


    # print(candidates)
    file.write(max(candidates, key=candidates.get) + "\n")
    # print(max(candidates, key=candidates.get))

file.close()