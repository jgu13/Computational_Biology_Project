import numpy as np
from ConsensusCellMask import segmentation_merge

def test_case1():
    seg_mapA = np.array(
        [[0, 1, 0, 0, 0],
         [1, 1, 1, 0, 0],
         [3, 1, 2, 2, 2],
         [3, 3, 2, 2, 2],
         [3, 3, 3, 2, 2],
         [3, 3, 0, 0, 0]])
    seg_mapB = np.array(
        [[0, 2, 0, 0, 0],
         [2, 2, 2, 0, 0],
         [3, 2, 1, 1, 1],
         [3, 3, 1, 1, 1],
         [3, 3, 3, 1, 1],
         [3, 3, 0, 0, 0]]
    )
    prob_mapA = np.array(
        [[0.10, 0.90, 0.05, 0.01, 0.00],
         [0.88, 0.03, 0.98, 0   , 0   ],
         [0.68, 0.90, 0.87, 0.97, 0.85],
         [0.66, 0.70, 0.85, 0.03, 0.80],
         [0.50, 0.03, 0.30, 0.90, 0.97],
         [0.50, 0.45, 0.40, 0   , 0   ]]
    )
    prob_mapB = np.array(
        [[0.00, 0.50, 0.00, 0.00, 0.00],
         [0.40, 0.03, 0.38, 0, 0],
         [0.88, 0.22, 0.35, 0.25, 0.33],
         [0.89, 0.88, 0.44, 0.05, 0.30],
         [0.90, 0.10, 0.87, 0.33, 0.44],
         [0.98, 0.99, 0, 0, 0]]
    )
    nucleus = np.array(
        [[0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0]]
    )
    seg_maps = [seg_mapA, seg_mapB]
    prob_maps = [prob_mapA, prob_mapB]
    return segmentation_merge(seg_maps, prob_maps, nucleus)

def test_case2():
    prob_matA_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\21RD_CD20_Orig_Probabilities.npy"
    prob_matA = np.load(prob_matA_file)

    prob_matB_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\21RD_CD3_Orig_Probabilities.npy"
    prob_matB = np.load(prob_matB_file)

    seg_mapA_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\Cell_Object_Image.npy"
    seg_matA = np.load(seg_mapA_file)

    seg_mapB_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\Cell_Object_Image.npy"
    seg_matB = np.load(seg_mapB_file)

    nuc_bin_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\NucleiMask.npy"
    nuc_bin = np.load(nuc_bin_file)

    seg_maps = [seg_matA, seg_matB]
    prob_maps = [prob_matA, prob_matB]
    segmentation_merge(seg_maps, prob_maps, nuc_bin)

def main():
    # seg_maps = list([np.array([[1,2,1],[2,1,2],[1,1,1]]), np.array([[1,2,1],[2,1,2],[1,1,1]])])
    # prob_maps = list([np.array([[0.90,0.87,0.99],[0.90,0.00,0.88],[0.90,0.89,0.99]]),np.array([[0.10,0.15,0.22],[0.20,0,0.22],[0.10,0.10,0.10]])])
    # print(seg_maps[0][0])
    # merged_mat = segmentation_merge(seg_maps, prob_maps)
    # print(merged_mat)
    try:
        # res_tuple = test_case1()
        # print(res_tuple[0])
        # print(res_tuple[1])
        # print(res_tuple[2])
        test_case1()
        test_case2()
    except Exception as error:
        print("Error occured in testcase1: " + str(error))

if __name__ == "__main__":
    main()