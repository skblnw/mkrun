#!/usr/bin/env python
__description__ = \
"""
Take a file containing frame number vs. distance, then identify frames that 
sample those distances at approximately sample_interval.  Optionally takes an
arbitrary number of template files.  The program searches the contents of these
files for a search string (by default XXX), replaces the search string with the
frame number, then writes out each file with a unique, frame-specific filename.
This last feature means that this script can be used to automatically generate
input files to run each umbrella sampling simulation. 

v. 0.1: fixed bug where the script would choke on distance_files that went by
        frame intervals besides steps of 1.
"""
__author__ = "Michael J. Harms [harmsm@gmail.com]"
__date__ = "120314"
__usage__ = "setupUmbrella.py distance_file sample_interval [template_file1 template_file2 ...]"
__version__ = "0.1"

import sys, os, re

def readDistanceFile(distance_file):
    """
    Read a distance file that correlates file number with center-of-mass 
    distance between pulling groups.  Format is:

    0 com0
    1 com1
    2 com2
    ...

    If this data was generated using distances.pl (Justin Lemkul), it could
    have multiple com calculations in the same file.  This function starts from
    the bottom of the file and only takes the last output.
    """

    # Read file
    f = open(distance_file,'r')
    lines = f.readlines()
    f.close()

    # Read the data from the bottom
    out_dict = {}
    for i in range(len(lines)-1,-1,-1):

        # Split on white-space; grab frame/distance
        columns = lines[i].split()
        key = int(columns[0])
        value = float(columns[1])

        # If we've already seen this key, we're done reading the unique data
        # from the file.
        if key in out_dict:
            break     
        else:
           out_dict[key] = value

    # Now put the values into a simple list, sorting to make sure that they are
    # in the correct order.
    keys = out_dict.keys()
    keys.sort()
    out = [(k,out_dict[k]) for k in keys]

    return out


def sampleDistances(distance_table,sample_interval):
    """
    Go through the distances list and sample frames at sample_interval.  
    Appropriate samples are identified by looking forward through the
    distances to find the one that is closest to current_distance +
    sample_interval.  
    """
   
    distances = [d[1] for d in distance_table]

    current_index = 0
    sampled_indexes = [current_index]
    while current_index < len(distances):
        
        target_distance = distances[current_index] - sample_interval
     
        # Walk through the rest of the distances list and find the distance 
        # that is closest to the target distance
        onward = [abs(target_distance-d) for d in distances[current_index:]]
        next_index = onward.index(min(onward)) + current_index
      
        # If we run out of distances to compare, our next_index will be the
        # current index.  This means we're done.  Otherwise, record that we
        # found a new index.
        if current_index == next_index:
            break
        else:
            sampled_indexes.append(next_index)
            current_index = next_index

    return sampled_indexes


def createOutputFile(template_file,frame_number,search_string="XXX"):
    """
    Look for instances of the search string in a template file and replace with
    the frame number in a new file.
    """

    out_file = "frame-%i_%s" % (frame_number,template_file)

    # Prompt the user before wiping out an existing file
    if os.path.exists(out_file):
        answer = raw_input("%s exists!  Overwrite (y|n)?" % out_file)
        answer = answer[0].lower()
        if answer != "y":
            return None
    
    # Read the contents of the template file
    f = open(template_file,'r')
    file_contents = f.read()
    f.close()

    # Write out the template file contents, replacing all instances of 
    # the search string with the frame number
    f = open(out_file,'w')
    f.write(re.sub(search_string,"%i" % frame_number,file_contents))
    f.close()    


def main(argv=None):
    """
    Parse command line, etc.
    """

    if argv == None:
        argv = sys.argv[1:]

    # Parse required command line arguments
    try:
        distance_file = argv[0]
        sample_interval = float(argv[1])
    except (IndexError,ValueError):
        err = "Incorrect command line arguments!\n\n%s\n\n" % __usage__
        raise IOError(err)

    # See if a template file has been specified
    try:
        template_files = argv[2:]
    except IndexError:
        template_files = []

    # Figure out which frames to use
    distance_table = readDistanceFile(distance_file)
    sampled_indexes = sampleDistances(distance_table,sample_interval)

    # If any template files were specified, use them to make frame-specific
    # output
    if len(template_files) != 0:
        print "Creating frame-specific output for files:"
        print "\n".join(template_files)
        for t in template_files:
            for i in sampled_indexes:
                frame = distance_table[i][0]
                createOutputFile(t,frame,search_string="XXX")

    # Print out summary of the frames we identified      
    out = ["%10s%10s%10s\n" % ("frame","dist","d_dist")]
    for i in range(len(sampled_indexes)):

        frame = distance_table[sampled_indexes[i]][0]
        dist = distance_table[sampled_indexes[i]][1]
        if i == 0:
            delta_dist = "%10s" % "NA" 
        else:
            prev_dist = distance_table[sampled_indexes[i-1]][1]
            delta_dist = "%10.3f" % (dist - prev_dist)

        out.append("%10i%10.3f%s\n" % (frame,dist,delta_dist))

    return out
   

# If called from the command line... 
if __name__ == "__main__":
    out = main() 
    print "".join(out)
