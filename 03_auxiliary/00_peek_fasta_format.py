import gzip, itertools

path = r"E:\华为\deepseek harness\gene_famaily_analysis\Araport11_pep_20250411.gz"
with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
    for i, line in enumerate(itertools.islice(fh, 12)):
        print(repr(line[:200]))
