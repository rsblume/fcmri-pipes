
## Functions
"""
# smaller are better, easier to undertsand and debug

 more than @ 5 inputs suggests bad data model, might want to 
  1. break into more pieces
  2. consider logic

For example lets look at run_correlations

1. too many inputs,and not clear what they are, if you do start with this make
sure to add explanation to docstring explaining what the inputs represent


Lets break into parts
"""
def command_line(command):
    """ run a shell command using subprocess"""
    proc = subprocess.Popen(command, stdout = subprocess.PIPE, shell=True)
    output, stderr = proc.communicate()
    retcode = proc.returncode
    if not retcode == 0:
        print stderr
        return None
    else:
        return output


def load_mask(roi2roidir, subid, nnodes, atlas, min_dist):
   """ load mask defined by
   roi2roidir : string that points to location of directory
   subid :  string representing subject id (eg 204)
   atlas: atlas used ('aal', 'dosenbach', 'dosenbach_all', 'HOall',
   'Greicius_func')
   min_dist : string defining/loading specific distance mask
   """
   try:
       maskfilename = os.path.join(roi2roidir, 
                                   'anatomdist_masks',
                                   '%s_%s_mindist$s.txt'%(subid, atlas))
       mask = np.load(maskfilename)
   except:
       print 'Warning, %s not found, using ones mask'%maskfilename
       mask = np.ones((nnodes, nnodes))
   return mask


def good_overlap(roi, masked_roi, thr= 0.25):
    """ given two rois (one masked) calculate the difference in volume
    if less than thr, return False, else return True
    """
    outputs = []
    for item in (roi, masked_roi):
        out = commandline(' '.join(['fslstats', roi, '-V'])
        outputs.append(float(out.split()[0]))
    proportion = outputs[1] / outputs[0] # maskedroi / roi
    if proportion >= thr:
        return True
    return False

