import numpy as np
import operator as op
import matplotlib.pyplot as plt
import itertools
import random as random

# class Pixel:
#     def __init__(self, i, j, prob, obj_id, marker):
#         self.i = i
#         self.j = j
#         self.obj_id = obj_id
#         self.prob = prob
#         self.marker = marker

class Cell:
    def __init__(self, marker, id, prob):
        self.marker = marker
        self.id = id
        self.prob = prob
        self.membrane = list() #list of coordinates
        self.nucleus = list() #list of coordinates

def dimention_check(nrow, ncol, rows, cols):
    return rows == nrow and cols == ncol

def average_prob(seg_map, prob_map):
    nrow = seg_map.shape[0]
    ncol = seg_map.shape[1]
    
    obj_prob_sum = [0]*(nrow*ncol)
    pix_count = [0]*(nrow*ncol)
    
    for i in range(nrow):
        for j in range(ncol):
            id = int(seg_map[i][j])
            if(id > 0):
                obj_prob_sum[id] += prob_map[i][j]
                pix_count[id] += 1
    
    for i in range(nrow):
        for j in range(ncol):
            id = int(seg_map[i][j])
            if(id > 0):
                prob_map[i][j]=obj_prob_sum[id] / pix_count[id]

def outline_cell(cell, res, seg, outline_res):
    nrow = len(res)
    ncol = len(res[0])
    #iterate thought the membrane pixel of this cell
    for p in cell.membrane:
        r = p[0]
        c = p[1]
        #if this pixel has res[pix[i][j]]==c.marker and its neighbor has a different id from itself
        if (res[r][c]==cell.marker):
            for i,j in zip([-1,1,0,0,-1,-1,1,1],[0,0,-1,1,-1,1,-1,1]):
                nr = r + i
                nc = c + j
                if (nr >= 0 and nr < nrow and nc >= 0 and nc < ncol and seg[nr][nc] != cell.id):
                    #then mark this pixel as outline
                    outline_res[r][c] = True

def colored(marker, prob, outline, num_markers):
    nrow, ncol = marker.shape[0], marker.shape[1]
    #initialize a color matrix
    img_mat = np.zeros((nrow, ncol,3),dtype = np.ubyte)
    #colors to be used
    colors = [(0,0,0)]*(num_markers+1)
    for i in range(num_markers):
        r = random.randint(0,255)
        g = random.randint(0,255)
        b = random.randint(0,255)
        colors[i+1] = (r,g,b)
    #color the matrix pixel by pixel
    for h in range(nrow):
        for w in range(ncol):
            if(outline[h][w]):
                pix_color = marker[h][w]
                img_mat[h][w] = [colors[pix_color][0]*100%256, colors[pix_color][1]*100%256, colors[pix_color][2]*100%256]
            elif(marker[h][w] == 0):
                img_mat[h][w] =[0,0,0]
            else:
                pix_color = marker[h][w]
                pix_prob = prob[h][w]
                img_mat[h][w] = [colors[pix_color][0], colors[pix_color][1], colors[pix_color][2]]
                img_mat[h][w] = img_mat[h][w]*pix_prob
    #flip the image
    img_mat_copy = img_mat.copy()
    #if nrow is odd, divide the last row index by half
    half = nrow/2 if(nrow%2 == 0) else (nrow-1)/2
    for y in range(nrow):
        for x in range(ncol):
            if y < half:
                flipped = int((y + abs(half - y) * 2) - 1)
                #print("Original: "+str(img_mat[y][x])+" Flipped:"+str(img_mat_copy[flipped][x]))
                img_mat[y][x] = img_mat_copy[flipped][x]
            elif y > half:
                flipped = int((y - abs(half - y) * 2) - 1)
                img_mat[y][x] = img_mat_copy[flipped][x]
    #render the image matrix to an image
    plt.imshow(img_mat, interpolation='nearest', origin = 'lower')
    plt.show()

def segmentation_merge(seg_maps, prob_maps, nuc_bin):

    # check dimentions
    nrow = seg_maps[1].shape[0]
    ncol = seg_maps[1].shape[1]
    
    for i in range(len(seg_maps)):
        if ( not dimention_check(nrow, ncol, seg_maps[i].shape[0], seg_maps[i].shape[1]) or not dimention_check(nrow, ncol, prob_maps[i].shape[0], prob_maps[i].shape[1])):
            raise Exception("Some matrix in the list has different dimentions!")
        
    #Init result matrix
    res = np.full((nrow, ncol), 0)
    prob_res = np.zeros((nrow,ncol))
    outline_res = np.full((nrow, ncol),False)
    #Init vector of pixels
    vCell = list()
    
    #average probability
    #store probability mean of each object in each pixel
    for (s,p,i) in zip(seg_maps, prob_maps, range(len(prob_maps))):
        #replace each prob_map with a copy
        prob_maps[i] = p.copy()
        average_prob(s,prob_maps[i])
    
    #for every img, add its pixels to the vector of pixels
    marker = 0
    for (s,p) in zip(seg_maps, prob_maps): 
        marker = marker + 1
        # a set of checked id
        id_checked = set()
        #set of cells of the same type
        cell_list = set()
        for i in range(nrow):
            for j in range(ncol):
                if(s[i][j] != 0):
                    #if this cell has not been checked
                    if(s[i][j] not in id_checked):
                        #we create a new cell
                        cell = Cell(marker,s[i][j],p[i][j])
                        # add this cell to the list of cell objects
                        cell_list.add(cell)
                        # add its id to the checked_id list
                        id_checked.add(cell.id)
                    else:
                        #find the first cell that has the current id
                        cell = next(c for c in cell_list if c.id == s[i][j])
                    #check if this pixel belongs to a nucleus
                    if(nuc_bin[i][j]):
                        #add this pixel to the list of pixels of nucleus this cell
                        cell.nucleus.append((i,j))
                        #otherwise, add this pixel to the list of pixels of membrane
                    else:
                        cell.membrane.append((i,j))
         #merge the set of cells into the list all cells
        for c in cell_list:
            vCell.append(c)

    #sort the pixels in reverse order
    vCell = sorted(vCell, key=op.attrgetter('prob', 'id'), reverse = True)
    
    #loop through cells
    k = 0
    while(k < len(vCell)):
        #check conflict for current object
        c = vCell[k]
        obj_id = int(c.id)
        nuc_taken = False
        #check if the nucleus of the cell can fit in the final image
        for n_pixel in c.nucleus:
            if res[n_pixel[0]][n_pixel[1]] != 0:
                nuc_taken = True
                break
        #skip this cell if the nucleus does not fit
        if(nuc_taken):
            k += 1
            continue
        #otherwise, we add this cell to the final image
        else:
            k += 1
            for n_pixel in c.nucleus:
                #intensity = p[n_pixel[0]][n_pixel[1]]
                res[n_pixel[0]][n_pixel[1]] = c.marker #mark nuclei pixel
                prob_res[n_pixel[0]][n_pixel[1]] = c.prob #mark pixel with cell prob mean
            for m_pixel in c.membrane:
                # if this pixel conflict with those in the final image, retain without replacement
                if res[m_pixel[0]][m_pixel[1]] == 0:
                    #intensity = p[m_pixel[0]][m_pixel[1]]
                    res[m_pixel[0]][m_pixel[1]] = c.marker #mark membrane pixel multiplied by intensity?
                    prob_res[m_pixel[0]][m_pixel[1]] = c.prob #mark pixel with cell prob mean
            #outline this cell
            outline_cell(c, res, seg_maps[c.marker-1], outline_res)
    #return the result matrix 
    return colored(res, prob_res, outline_res, len(seg_maps))

