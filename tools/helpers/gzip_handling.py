import sys


def gzip_test(in1):
    try:
        with open(in1) as f:
            f.readline()
        compressed = False
    except UnicodeDecodeError:
        compressed = True
    except IsADirectoryError:
        return None
    return compressed


if __name__ == '__main__':
    gzip_test[sys.argv[1]]