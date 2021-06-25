import re
import sys
import traceback

import cv2 as cv
import numpy as np
import wx
import ConsensusCellMask as CCM


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
        [[False, False, False, False, False],
         [False, True, False, False, False],
         [False, False, False, False, False],
         [False, False, False, True, False],
         [False, True, False, False, False],
         [False, False, False, False, False]]
    )
    seg_maps = [seg_mapA, seg_mapB]
    prob_maps = [prob_mapA, prob_mapB]
    return CCM.segmentation_merge(seg_maps, prob_maps, nucleus)

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
    CCM.segmentation_merge(seg_maps, prob_maps, nuc_bin)

def test_case3():
    prob_matA_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\21RD_CD20_Orig_Probabilities.npy"
    prob_matA = np.load(prob_matA_file)

    prob_matB_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\21RD_CD3_Orig_Probabilities.npy"
    prob_matB = np.load(prob_matB_file)

    prob_matC_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD4\\21RD_CD4_Orig_Probabilities.npy"
    prob_matC = np.load(prob_matC_file)

    seg_mapA_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\Cell_Object_Image.npy"
    seg_matA = np.load(seg_mapA_file)

    seg_mapB_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\Cell_Object_Image.npy"
    seg_matB = np.load(seg_mapB_file)

    seg_mapC_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD4\\Cell_Object_Image.npy"
    seg_matC = np.load(seg_mapC_file)

    nuc_bin_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\NucleiMask.npy"
    nuc_bin = np.load(nuc_bin_file)

    seg_maps = [seg_matA, seg_matB, seg_matC]
    prob_maps = [prob_matA, prob_matB, prob_matC]
    CCM.segmentation_merge(seg_maps, prob_maps, nuc_bin)

def test_case4():
    prob_matA_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\21RD_CD20_Orig_Probabilities.npy"
    prob_matA = np.load(prob_matA_file)

    prob_matB_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\21RD_CD3_Orig_Probabilities.npy"
    prob_matB = np.load(prob_matB_file)

    prob_matC_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD4\\21RD_CD4_Orig_Probabilities.npy"
    prob_matC = np.load(prob_matC_file)

    prob_matD_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD8a\\21RD_CD8a_Orig_Probabilities.npy"
    prob_matD = np.load(prob_matD_file)

    seg_mapA_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\Cell_Object_Image.npy"
    seg_matA = np.load(seg_mapA_file)

    seg_mapB_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\Cell_Object_Image.npy"
    seg_matB = np.load(seg_mapB_file)

    seg_mapC_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD4\\Cell_Object_Image.npy"
    seg_matC = np.load(seg_mapC_file)

    seg_mapD_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD8a\\Cell_Object_Image.npy"
    seg_matD = np.load(seg_mapD_file)

    nuc_bin_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\NucleiMask.npy"
    nuc_bin = np.load(nuc_bin_file)

    seg_maps = [seg_matA, seg_matB, seg_matC, seg_matD]
    prob_maps = [prob_matA, prob_matB, prob_matC, prob_matD]
    return CCM.segmentation_merge(seg_maps, prob_maps, nuc_bin)

def get_target_index(target):
    p = re.compile(target)
    for i, f in zip(range(len(files)), files):
        target = p.search(f)
        if (target != None):
            return i
    return None

def get_targets(files):
    targets = []
    p = re.compile('CD[0-9a-zA-Z]+')
    for f in files:
        target = p.search(f)
        if target != None:
            targets.append(target.group())
    return targets

def overlay_test1(res, outline, prob_res, dest=None):
    app = wx.App()
    fileDialog = wx.FileDialog(None, 'Open', wildcard="(*.tiff)|*.tiff|(*.png)|*.png|(*.jpg)|*.jpg",
                               style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    if fileDialog.ShowModal() == wx.ID_OK:
        file_path = fileDialog.GetPath()
    else:
        print("Image not loaded. Abort overlay_test1.")
        return
    img = cv.imread(file_path)
    target = input("Please enter the marker signal to be outlined: e.g. CD20 ")
    index = get_target_index(target)
    if index != None:
        CCM.overlay(img,outline,res,index+1, default_dest=dest) #marker is index+1
    else:
        raise("Target not found.")

def get_standard_prob_mean_test():
    res, outline_res, prob_res, res_cells_list = test_case4()
    CCM.get_standard_prob_mean(res_cells_list)

def Open(type):
    #create an object of application class
    app = wx.App()
    fileDialog = wx.FileDialog(None, str('Open'+type), wildcard=str("*" + type + "*.npy"),
                               style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    if fileDialog.ShowModal() == wx.ID_OK:
        path = fileDialog.GetPath()
    else:
        path = None
    fileDialog.Destroy()
    return path

def prompt_for_images():
    stop = False
    while not stop:
        #open probability npy file
        target_file = Open("Probabilities")
        if(target_file == None):
            raise("No file selected.")
        target_prob_mat = np.load(target_file)
        p = re.compile(r"\d{2}\w{2}_\w{2}[0-9A-Za-z]*(_Orig)?_Probabilities.npy")
        f = p.search(target_file)
        if(f==None):
            raise("File name does not match.")
        else:
            files.append(f.group()) #store matched strings in files
        #open segmentation npy file
        target_file = Open("Cell_Object_Image")
        target_seg_mat = np.load(target_file)

        seg_maps.append(target_seg_mat)
        prob_maps.append(target_prob_mat)
        ans = input("Do you want to enter more targets? Enter y for yes, n for no ")
        while True:
            ans = ans.lower()
            if (ans == "n"):
                stop = True
                break
            elif (ans == "y"):
                break
            else:
                ans = input("You did not enter the right answer. Enter Y for yes, N for no ")
    # open nuclei bin file
    target_file = Open("NucleiMask")
    nuc_bin = np.load(target_file)
    #extract sample name
    end = 0
    for i in range(len(files[0])):
        if files[0][i] == '_':
            end = i #excluding '_'
            break
    sample = files[0][0:end]
    return seg_maps, prob_maps, nuc_bin, files, sample

def show_percent_coverage(prob_maps, res, sample):
    #TODO:
    # For each merged sample, plot percent coverage
    targets = get_targets(files)
    CCM.plot_percent_coverage(prob_maps, res, targets, sample)

def main():
    try:
        # res_tuple = test_case1()
        # print(res_tuple[0])
        # print(res_tuple[1])
        # print(res_tuple[2])
        # test_case1()
        # test_case2()
        # test_case3()
        # test_case4()
        # get_standard_prob_mean_test()

        global seg_maps
        global prob_maps
        global nuc_bin
        global files

        seg_maps = []
        prob_maps = []
        files = []
        res_cell_lists = {}
        done = False

        # nuc_bin_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\NucleiMask.npy"
        # nuc_bin = np.load(nuc_bin_file)

        #TODO:
        # while not done, keep asking user for merging sample cell masks
        while not done:
            if(len(seg_maps)>0):
                ans = input("Do you want to keep all previously selected images? Enter Y for yes, N for no: ")
                while ans.lower() != 'y' and ans.lower() != 'n':
                    ans = input("Your answer is invalid. Enter Y for yes, N for no: ")
                if(ans == "n"):
                    seg_maps = []
                    prob_maps = []
                    files = []
                    seg_maps, prob_maps, nuc_bin, files, sample = prompt_for_images()
            else:
                seg_maps, prob_maps, nuc_bin, files, sample = prompt_for_images()

            ans = input("Would you like to save the merged cell masks? Enter Y for yes, N for no: ")
            while ans.lower()!='y' and ans.lower()!='n':
                ans = input("Your ans is invalid. Enter Y for yes, N for no: ")
            if ans.lower() == 'y':
                dest = input("Enter a location followed by a file name.\n(For example: foo\\bar\\merged_mask.tiff) ")
                res, outline_res, prob_res, res_cells_list = CCM.segmentation_merge(seg_maps, prob_maps, nuc_bin, dest)
            else:
                res, outline_res, prob_res, res_cells_list = CCM.segmentation_merge(seg_maps, prob_maps, nuc_bin)
            res_cell_lists[sample] = (res_cells_list, files)
            show_percent_coverage(prob_maps, res, sample)
            print("{} cells identified in {}".format(len(res_cells_list), sample))

            while True:
                ans = input("Would you like to overlay outline? Enter Y for yes, N for no: ")
                ans = ans.lower()
                if(ans == "y"):
                    ans=input("Would you like to save the result? Enter Y for yes, N for no: ")
                    while ans.lower() != 'y' and ans.lower() != 'n':
                        ans = input("Your answer is invalid. Enter Y for yes, N for no.")
                    if (ans.lower() == 'y'):
                        dest = input("Please enter a location followed by the name you want to save the image. \n For instance: foo\\bar\\CD20outline.png")
                    elif(ans.lower() == 'n'):
                        dest = None

                    overlay_test1(res, outline_res, prob_res, dest)
                elif(ans == "n"):
                    ans = input("Do you want to continue merging cell masks? Enter Y for yes, N for no:")
                    while ans.lower() != 'y' and ans.lower() != 'n':
                        ans = input("Your answer is invalid. Enter Y for yes, N for no.")
                    if (ans.lower() == 'y'):
                        break
                    elif(ans.lower() == 'n'):
                        done=True
                        break

        #TODO:
        # For all merged samples, plot cell area and error rates distribition
        file_lists = []

        for key in list(res_cell_lists.keys()):
            file_lists.append(list(res_cell_lists[key][1]))

        #assume users selected same markers for all samples
        markers = get_targets(file_lists[0])

        CCM.plot_cells_area(list(res_cell_lists[key][0] for key in list(res_cell_lists.keys())),
                            list(res_cell_lists.keys())
                            )
        CCM.plot_error_rates(list(res_cell_lists[key][0] for key in list(res_cell_lists.keys())),
                             markers,
                             list(res_cell_lists.keys())
                             )
    except Exception as error:
        print("Error occured in main: " + str(error))
        traceback.print_exc(file=sys.stdout)

if __name__ == "__main__":
    main()