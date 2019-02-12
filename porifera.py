import sys


def main():
    in1 = sys.argv[1]
    max_len = porifera(in1)


def porifera(in1):
    with open(in1) as f:
        max_len, y = 0, 0
        for line in f:
            y += 1
            if y == 2:
                len_test = len(line.rstrip())
                if len_test > max_len:
                    max_len = len_test
    return max_len


if __name__ == "__main__":
    main()
