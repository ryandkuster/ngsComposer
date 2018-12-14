def main():
    barcodes_matrix = []
    with open('key.txt') as f:
        R2_barcodes = R2_barcodes_maker([], f)
        barcodes_matrix = array_maker(R2_barcodes)
        barcodes_matrix, R1_barcodes = barcodes_matrix_maker([], f, barcodes_matrix)
        print(barcodes_matrix)
        print(R1_barcodes)
        print(R2_barcodes)


def barcodes_matrix_maker(R1_barcodes, f, barcodes_matrix):
    for i, line in enumerate(f):
        for j, item in enumerate(line.split()):
            if j == 0:
                R1_barcodes.append(item)
            if j > 0:
                barcodes_matrix[j - 1].append(item)    
    return barcodes_matrix, R1_barcodes


def R2_barcodes_maker(R2_barcodes, f):
    line = f.readline()
    for item in line.split():
        R2_barcodes.append(item)
    return R2_barcodes

        
def array_maker(R2_barcodes):
    barcodes_matrix = [[0] * 1 for i in range(len(R2_barcodes))]
    for i, x in enumerate(R2_barcodes):
        barcodes_matrix[i][0] = x
    return barcodes_matrix

if __name__ == "__main__":
    main()