import numpy as np
import operator as op
import matplotlib.pyplot as plt
import random as random

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

def colored(marker, prob, outline, num_markers, default_dest):
    nrow, ncol = marker.shape[0], marker.shape[1]
    #initialize a color matrix
    img_mat = np.zeros((nrow, ncol,3),dtype = np.ubyte)
    #colors to be used
    colors = [(0,0,0)]*(num_markers+1)
    for i in range(num_markers):
        random.seed(i*4)
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
    if(default_dest != None):
        plt.imsave(default_dest, img_mat)
    #clear before drawing new image
    fig = plt.figure()
    fig.tight_layout()
    # render the image matrix to an image
    plt.imshow(img_mat)
    plt.show()

def binary_search(list, key):
    length = len(list)
    if(length < 0 ):
        raise("Error in binary search: Length less than 0.")
    if(length == 1):
        return list[0]
    low, high = 0, length-1
    while low <= high:
        mid = int((high - low) / 2) + low
        if(list[mid].id == key):
            return list[mid]
        elif(list[mid].id > key):
            high = mid - 1
        else: #list[mid].id < key
            low = mid + 1
    return None

def segmentation_merge(seg_maps, prob_maps, nuc_bin, default_dest=None):

    # check dimentions
    nrow = seg_maps[0].shape[0]
    ncol = seg_maps[0].shape[1]
    
    for i in range(len(seg_maps)):
        if ( not dimention_check(nrow, ncol, seg_maps[i].shape[0], seg_maps[i].shape[1]) or not dimention_check(nrow, ncol, prob_maps[i].shape[0], prob_maps[i].shape[1])):
            raise Exception("Some matrix in the list has different dimentions!")
        
    #Initiation
    res = np.full((nrow, ncol), 0)
    prob_res = np.zeros((nrow,ncol))
    outline_res = np.full((nrow, ncol),False)
    res_cells_list = set()
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
        cell_list = list()
        for i in range(nrow):
            for j in range(ncol):
                if(s[i][j] != 0):
                    #if this cell has not been checked
                    if(s[i][j] not in id_checked):
                        #we create a new cell
                        cell = Cell(marker,s[i][j],p[i][j])
                        # add this cell to the list of cell objects
                        cell_list.append(cell)
                        # add its id to the checked_id list
                        id_checked.add(cell.id)
                    else:
                        #find the first cell that has the current id
                        #cell = next(c for c in cell_list if c.id == s[i][j]) #linear search taking forever
                        #sort cell list by increasing id
                        cell_list = sorted(cell_list, key=op.attrgetter('id'))
                        #binary_search
                        cell = binary_search(cell_list, s[i][j])
                        if (cell == None):
                            raise("No cell with id: {} found.".format(s[i][j]))
                    #check if this pixel belongs to a nucleus
                    if(nuc_bin[i][j]):
                        #add this pixel to the list of pixels of nucleus this cell
                        cell.nucleus.append((i,j))
                        #otherwise, add this pixel to the list of pixels of membrane
                    else:
                        cell.membrane.append((i,j))
         #merge the set of cells into the list all cells
        vCell.extend(cell_list)

    #sort the pixels in reverse order
    vCell = sorted(vCell, key=op.attrgetter('prob'), reverse = True)
    
    #loop through cells
    k = 0
    count = 0
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
                    res[m_pixel[0]][m_pixel[1]] = c.marker
                    prob_res[m_pixel[0]][m_pixel[1]] = c.prob #mark pixel with cell prob mean
            #outline this cell
            outline_cell(c, res, seg_maps[c.marker-1], outline_res)
            #re-number and save this cell into the res_cells_list. Since one cell will not be considered more than once, it is safe to re-number
            c.id = count+1
            res_cells_list.add(c)
    #return the result matrix 
    colored(res, prob_res, outline_res, len(seg_maps), default_dest)
    return res, outline_res, prob_res, res_cells_list

#TODO:
# overlay cell outline onto colored marker signal image
def overlay(signal_img, outline, marker, target, default_dest = None):
    #check if all signal imgs have the same dimentions as outline
    nrow, ncol = outline.shape[0], outline.shape[1]
    #check for dimension
    if(not dimention_check(nrow, ncol, signal_img.shape[0], signal_img.shape[1])):
        raise("Wrong dimention: Image is {}x{}. target is {}x{}.".format(signal_img.shape[0],signal_img.shape[1],nrow, ncol))

    outlined = np.copy(signal_img)
    #if the pixel is marker and outline is true, marks it as white
    for i in range(nrow):
        for j in range(ncol):
            if (marker[i][j]==target and outline[i][j]):
                outlined[i][j] = [255,255,255]
    plt.figure()
    plt.imshow(outlined,interpolation='none')
    if(default_dest!=None):
        plt.imsave(default_dest, outlined)

def get_percent_covarage(prob_map, marker, target):
    #check dimentions!
    nrow, ncol = marker.shape[0], marker.shape[1]
    # check for dimension
    if (not dimention_check(nrow, ncol, prob_map.shape[0], prob_map.shape[1])):
        raise (
            "Wrong dimention: Image is {}x{}. target is {}x{}.".format(prob_map.shape[0], prob_map.shape[1], nrow, ncol))

    sum, high_sum, coverage, high_coverage= 0, 0, 0, 0
    #TODO:
    # iterate thought the prob_map to count pixels having probability > 0,
    # if this pixel is also the target, count as coverage
    for i in range(nrow):
        for j in range(ncol):
            if prob_map[i][j] > 0:
                sum += 1
                if prob_map[i][j] > 0.5:
                     high_sum += 1
                     if marker[i][j] == target:
                        high_coverage += 1
                if marker[i][j] == target:
                    coverage += 1

    #return the coverage/sum
    return (coverage/sum)*100, (high_coverage/high_sum)*100

# This function shows the signal coverage of each marker for each sample
def plot_percent_coverage(prob_maps, marker, targets, sample):
    #TODO:
    # for each marker, get its percent coverage
    # store it in percent coverage vector and plot with X = marker names, Y = percent coverage
    percentages = {}
    for prob_map, target, target_index in zip(prob_maps, targets, range(len(targets))):
        key = target
        low_val, high_val = get_percent_covarage(prob_map, marker, target_index+1)
        percentages[key] = (low_val, high_val)
    #clear canvas before plotting
    fig, ax = plt.subplots()
    data = []
    data.append([i[0] for i in list(percentages.values())])
    data.append([i[1] for i in list(percentages.values())])
    X = np.arange(len(targets)) #bar positions
    rec1 = ax.bar(X+0.00, height=data[0], color='b', width=0.35)
    rec2 = ax.bar(X+0.35, height=data[1], color='g', width=0.35)
    ax.set_title("{} signal percentage coverage".format(sample))
    ax.set_xticks(X)
    ax.set_xticklabels(targets)
    ax.set_xlabel("Markers")
    ax.set_ylabel("Percent coverage")
    ax.legend(labels=['Total signal', 'High signal'])
    ax.bar_label(rec1, padding=3)
    ax.bar_label(rec2, padding=3)
    fig.tight_layout()

# This function gets all cell areas
def get_sample_cells_area(res_cells_list):
    #TODO:
    # return list of cell areas
    data = []
    for c in res_cells_list:
        data.append(len(c.membrane) + len(c.nucleus))
    return data

#This function plot cell areas of each sample as box plot
def plot_cells_area(res_cell_lists, samples):
    if len(res_cell_lists) != len(samples):
        raise("Length of cell lists and length of samples do not match!")

    data_list = []
    for i in range(len(samples)):
        sample_cell_area = get_sample_cells_area(res_cell_lists[i])
        data_list.append(sample_cell_area)

    fig = plt.figure(figsize=(6,7))
    ax = fig.add_subplot(111)

    bp = plt.boxplot(data_list, vert=True, labels=samples, meanline=True)

# This function is used once before getting error rates to get the lower bound of probability mean
def get_standard_prob_mean(test_res_cells_list):
    #TODO:
    # plot probability mean of all cells
    X = np.array([c.prob for c in test_res_cells_list])
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.hist(X)
    plt.title("Cell probability mean distribution")

# This function get the error rates for each marker of a sample
def get_sample_error_rates(res_cells_list, num_markers):
    #TODO:
    # On average, how many cells out of all cells have probability mean less than 25% for each cell type
    # For each cell, get its cell.prob
    # if the prob mean is less than 25%,then add it to num_cells and record its

    sample_error_rates = [0]*(num_markers+1)
    for c in res_cells_list:
        if c.prob < 0.25:
            sample_error_rates[c.marker-1] += 1
    for r in range(len(sample_error_rates)):
        sample_error_rates[r] = (sample_error_rates[r] / len(res_cells_list))*100

    #append total error rate at the end
    sample_error_rates[len(sample_error_rates)-1] = sum(sample_error_rates)
    return sample_error_rates

def plot_error_rates(res_cell_lists, markers, samples):
    #TODO:
    # tabulate sample_error_rates with rows as samples, cols as markers

    all_sample_error_rates = []
    for list in res_cell_lists:
        all_sample_error_rates.append(get_sample_error_rates(list, len(markers)))

    all_sample_error_rates = [['%.2f' % j for j in i] for i in all_sample_error_rates]
    colors = plt.cm.BuPu(np.linspace(0,0.5,len(samples)))
    markers.append('Total')

    fig, ax = plt.subplots()
    table = ax.table(cellText=all_sample_error_rates,
                     rowLabels=samples,
                     colLabels=markers,
                     cellLoc='center',
                     loc='upper left'
                     )

    table.set_fontsize(15)
    table.scale(1,2)
    ax.axis('off')
    plt.title('Sample error rates')
    plt.show()



