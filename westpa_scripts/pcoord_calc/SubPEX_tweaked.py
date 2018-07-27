#!/usr/bin/python

# import dependencies
import multiprocessing
import numpy as np
import scipy as sp
import scoria

def get_settings(*path):
    # Getter Function for in-file Paramerters.
    # Retruns [resolution, dimension, filepath,ligands]

    # Open File
    loc = "SubPEX_settings"                                             # Default Settings Location
    if(len(path) > 0):                                                  # If Non-Defaults
        loc = path[0]                                                   # Set Them
    subpex_settings = open(loc)                                         # Open File
                                                                        #
    settings = [None] * 7                                               # No Defaults 
                                                                        #
    for line in subpex_settings:                                        #
        fields = line.strip().split()                                   # Splits Columns
                                                                        #
        if(fields[0].lower() == "resolution"):                          # Set Resolution
            settings[0] = float(fields[1])                              #
        elif(fields[0].lower() == "dimension"):                         # Set Dimension
            settings[1] = float(fields[1])                              #
        elif(fields[0].lower() == "filepath"):                          # Set File Path
            settings[2] = fields[1]                                     #
        elif(fields[0].lower() == "center"):                            # Set Center Coordinates
            coords = fields[1].split(",")                               # Coordinates are Comma Delimited  
            settings[3] = float(coords[0])                              #
            settings[4] = float(coords[1])                              #
            settings[5] = float(coords[2])                              #
        elif(fields[0].lower() == "chain"):                             # Set Chain
            settings[6] = fields[1]                                     #
            global chain_id                                             #
            chain_id = fields[1].upper()                                #
        else:
            print "Error: Unknown Field in SubPEX Settings"
    
    return settings

def get_ligand_info(chain_id, ligand_name):
	ligand = mol.select_atoms({'chainid':chain_id,'resname':ligand_name})
	ligand_com = np.around(mol.get_center_of_mass(ligand),3)
	return ligand_com

# get_snapped_center is replacing the old CenterPoint class

def get_snapped_center(center, resolution):
	snapped_x = np.around(center[0]/resolution)*resolution
	snapped_y = np.around(center[1]/resolution)*resolution
	snapped_z = np.around(center[2]/resolution)*resolution
	snapped_center = np.around(np.array([snapped_x, snapped_y, snapped_z],float),3)
	return snapped_center

# field_of_points is replacing the old FieldOfPoints class

def field_of_points():
	resolution = inputs[0]
	dimension = inputs[1]
	xs = np.round(np.arange(center[0] - dimension/2, center[0] + dimension/2 + resolution, resolution),3)
	ys = np.round(np.arange(center[1] - dimension/2, center[1] + dimension/2 + resolution, resolution),3)
	zs = np.round(np.arange(center[2] - dimension/2, center[2] + dimension/2 + resolution, resolution),3)
	field_of_points = np.vstack(np.meshgrid(xs,ys,zs))
	field_of_points = field_of_points.reshape(3,-1).T
	return field_of_points

def get_all_alphas():
	alphas = mol.select_atoms({'name':['CA']})
	alphas = mol.get_molecule_from_selection(alphas)
	return alphas

# keep_pocket_alphas and keep_pocket_atoms are replacing the old KeepUsefulAlphas class

def keep_pocket_alphas():
	minimum = np.amin(reference_fop, axis=0)
	maximum = np.amax(reference_fop, axis=0)
	kept_pocket_alphas = []
	for frame in all_alpha_trajectory:
		keep = []
		for point in frame:
			if minimum[0] < point[0] < maximum[0]:
				if minimum[1] < point[1] < maximum[1]:
					if minimum[2] < point[2] < maximum[2]:
						keep.append(point)
		kept_pocket_alphas.append(np.array(keep))
	return np.array(kept_pocket_alphas)

def keep_pocket_atoms():
	minimum = np.amin(reference_fop, axis=0)
	maximum = np.amax(reference_fop, axis=0)
	kept_pocket_atoms = []
	atoms = mol.select_atoms({"chainid":[chain_id],"record_name":["ATOM  "]})
	atoms = mol.get_molecule_from_selection(atoms)
	atom_trajectory = atoms.get_trajectory_coordinates()
	for frame in atom_trajectory:
		keep = []
		for point in frame:
			if minimum[0] < point[0] < maximum[0]:
				if minimum[1] < point[1] < maximum[1]:
					if minimum[2] < point[2] < maximum[2]:
						keep.append(point)
		kept_pocket_atoms.append(np.array(keep))
	return np.array(kept_pocket_atoms)

def eliminate_points_that_clash_w_protein(point_field):
	for index in r:
		for item in index:
			point_field[item] = None

# the following several functions are courtesy of Patrick Ropp

def multithreading(inputs, num_processors, task_name):
    """Initialize this object.
    Args:
        inputs ([data]): A list of data. Each datum contains the details to
            run a single job on a single processor.
        num_processors (int): The number of processors to use.
        task_name (class): The class that governs what to do for each
            job on each processor.
    
    """
    results = []
    # If there are no inputs, just return an empty list.
    if len(inputs) == 0:
        return results
    num_processors = count_processors(len(inputs), num_processors)
    tasks = []
    for item in inputs:
        if not isinstance(item, tuple):
            item = (item,)
        tasks.append((task_name, item))
    results = start_processes(tasks, num_processors)
    return results

def flatten_list(tier_list):
    """
    Given a list of lists, this returns a flat list of all items.
    :params list tier_list: A 2D list.
    :returns: A flat list of all items.
    """
    flat_list = [item for sublist in tier_list for item in sublist]
    return flat_list


def strip_none(none_list):
    """
    Given a list that might contain None items, this returns a list with no
    None items.
    :params list none_list: A list that may contain None items.
    :returns: A list stripped of None items.
    """
    results = [x for x in none_list if x != None]
    return results
#
# Worker function
#

def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)

def count_processors(num_inputs, num_processors):
    """
    Checks processors available and returns a safe number of them to
    utilize.
    :param int num_inputs: The number of inputs.
    :param int num_processors: The number of desired processors.
    :returns: The number of processors to use.
    """
    # first, if num_processors <= 0, determine the number of processors to
    # use programatically
    if num_processors <= 0:
        num_processors = multiprocessing.cpu_count()
    # reduce the number of processors if too many have been specified
    if num_inputs < num_processors:
        num_processors = num_inputs
    return num_processors

def start_processes(inputs, num_processors):
    """
    Creates a queue of inputs and outputs
    """
    # Create queues
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()
    # Submit tasks
    for item in inputs:
        task_queue.put(item)
    # Start worker processes
    for i in range(num_processors):
        multiprocessing.Process(target=worker, args=(task_queue, done_queue)).start()
    # Get and print results
    results = []
    for i in range(len(inputs)):
        results.append(done_queue.get())
    # Tell child processes to stop
    for i in range(num_processors):
        task_queue.put('STOP')
    return results

############################## Start Analysis ##################################
# get inputs from user
#inputs = get_inputs()
settings = get_settings()
inputs = settings[0:3]


# initialize scoria Molecule by reading in the file from the inputs
mol = scoria.Molecule(inputs[-1])

# get the binding pocket center 
center_unsnapped = np.array(settings[3:6], float)
center = get_snapped_center(center_unsnapped, inputs[0])

# create a static field_of_points to be referenced but never modified
reference_fop = field_of_points()

# get all the alpha carbon indices
all_alphas = get_all_alphas()

# get the multi-frame coordinates of all the alpha carbons using aforementioned
# indices
all_alpha_trajectory = all_alphas.get_trajectory_coordinates()
print 'There are {} frames, each containing {} alpha carbons'.format(len(all_alpha_trajectory),len(all_alpha_trajectory[0]))

# only keep the alpha carbons coordinates in each frame that are in range of
# the binding pocket dimension parameter obtained from the user
kept_alphas = keep_pocket_alphas()

# summary outputs
for frame in kept_alphas:
	print "This frame has {} alpha carbons in the binding pocket".format(len(frame))

# get the coordinates for only those atoms whose coordinates in each frame are
# in range of the binding pocket dimension parameter obtained from the user
kept_atoms = keep_pocket_atoms()

# summary outputs
for frame in kept_atoms:
	print "This frame has {} atoms in the binding pocket".format(len(frame))

# now use the data to do the pocket analysis
for frame_index in xrange(len(kept_atoms)):                         # for each frame present grab
	frame_atoms = kept_atoms[frame_index]                           # coords for all atoms/alphas
	frame_alphas = kept_alphas[frame_index]                         # kept for this frame
	local_fop = reference_fop                                       # local dynamic copy of fop
	tree_atoms = sp.spatial.cKDTree(frame_atoms)                    # distances btw atoms
	tree_fop = sp.spatial.cKDTree(local_fop)                        # distances in field of points
	r = tree_atoms.query_ball_tree(tree_fop, 1.5, p=2, eps=0)       # points too close to atoms
	eliminate_points_that_clash_w_protein(local_fop)                # eliminate too close points
	hull = sp.spatial.ConvexHull(frame_alphas)                      # create reference hull

	def hull_checker_inputs():                                      # function to set up for hull checker
		hull_inputs = []                                        # because NAN points don't need to be
		for index in xrange(len(local_fop)):                    # checked... we already know they aren't
			if not np.isnan(index):                         # in the hull
				hull_inputs.append(index)
			else:
				pass
		return hull_inputs

	to_be_checked = hull_checker_inputs()                           # run the setup function

	def hull_checker(index):                                        # function to check whether a point
		point = reference_fop[index]                            # is within the convex hull made by
		if not np.isnan(point[0]):                              # kept alphas.  Essentially if adding
			points = list(frame_alphas)                     # that point to alphas changes the 
			points.append(point)                            # hull's vertices (see scipy.spatial
			points = np.array(points)                       # ConvexHull documentation), then it's 
			new_hull = sp.spatial.ConvexHull(points)        # outside the hull
			if list(hull.vertices) == list(new_hull.vertices):
				return index
			else:
				return None

        def chunk_data_for_multiprocess():               
                shape = []                                              # list to store indices of points in hull;
                if len(to_be_checked) <= 32766:
                        shape = multithreading(to_be_checked, 4, hull_checker)
                else:
                        chunk_indices = []
                        start = 0
                        end = 32767
                        while len(to_be_checked) > end - 1:
                                chunk_indices.append([start,end])
                                start += 32767
                                end += 32767
                        chunk_indices.append([start, len(to_be_checked)])
                        for i in chunk_indices:
                                in_hull = multithreading(a[i[0]:i[1]], 4, hull_checker)
                                for j in in_hull:
                                        shape.append(j)
                return shape

	pocket_shape = chunk_data_for_multiprocess()                    # if something other than None gets
        final_shape = []
        for index in pocket_shape:
                if index != None:
                        final_shape.append(reference_fop[index])
                else:
                        pass

	# summary outputs
	print len(final_shape)                                        
	print "Pocket volume: {}".format(len(final_shape)*(inputs[0]**3))

	# file output
	out_path = inputs[2][:-4] + '_frame{}'.format(frame_index) + '.xyz'
	g = open(out_path, 'w')
	g.write(str(len(final_shape)) + "\n")
	g.write("Generated in SubPEx.\tRes: {}\tDim: {} Ang\tVolume: {} Ang^3\n".format(inputs[0],inputs[1],len(final_shape)*(inputs[0]**3)))
	for point in final_shape:
		line = "CA\t{}\t{}\t{}\n".format(point[0],point[1],point[2])
		g.write(line)
	# summary output
	print "Done with frame {} of {}".format(frame_index+1, len(kept_atoms))
	
with open('pvol.txt', 'a+') as file:
    file.write("\n")
    file.write( "{}".format(len(final_shape)*(inputs[0]**3)))
